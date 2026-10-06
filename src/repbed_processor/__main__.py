import os
import click
import datetime
import pandas as pd
from pathlib import Path
from collections import Counter
from natsort import natsort_keygen
from repbed_processor.utils.logger import Logger


def filter_by_primers(df, primers, logger, keep_false):
    """
    filter to only the provided primers
    :param DataFrame df: dataframe of data for filtering
    :param primers: primers to filter by
    :param logger: logger object
    :param bool keep_false: whether to keep alignments with false in the primer suffix
    :return: filtered DataFrame
    """
    # filter to only relevant alignments (i.e. those whose qnames contain the specified primers)
    logger.info(f"Filtering to reads that end with the primer(s): {', '.join(primers)}")

    # initialise the list for dfs
    filtered_by_primer = []

    # get just the primer bit, to ensure no accidental matches of primers within qnames
    try:
        # qnames should be separated by .
        df["primer_suffix"] = (df["qname"].str.rsplit(".", n=1, expand=True))[1]
    except IndexError:
        try:
            # have seen some separated by /
            df["primer_suffix"] = (df["qname"].str.rsplit("/", n=1, expand=True))[1]
        except IndexError:
            logger.error(f"qnames in unexpected format. Could not check primers in qnames.")

    # for each primer, collect the data and add it to the list
    for primer_string in primers:
        primerfilterdf = df[df["primer_suffix"].str.contains(primer_string)]
        if not keep_false:
            logger.info(f"Discarding events with false primers for {primer_string}")
            primerfilterdf = primerfilterdf[~primerfilterdf["primer_suffix"].str.contains("false", case=False)]
        filtered_by_primer.append(primerfilterdf)

    # concatenate all the collected dfs
    primerfilterdf = pd.concat(filtered_by_primer)

    # where more than one primer is present in the read name, there can be repeated rows. remove them.
    final_version = primerfilterdf.drop_duplicates(keep='first')

    logger.info(f"Filtering done.")

    return final_version

def split_flanks_and_inserts(clusterdata, complex_limit, logger, primer, filter_inferred):
    """
    label the first and last alignments in each read (by alignment position on read)
    :param DataFrame clusterdata: df of data for clustering
    :param int complex_limit: number of complex alignments
    :param Logger logger: the logger object
    :param list primer: name of the primer(s)
    :return: data df, and lists of relevant captured data
    """
    # initialise lists of information to collect
    aln_lens = []
    flanking_aln_sizes = []
    insertion_aln_sizes = []
    simple_insert = []
    complex_insert = []
    chromflank = []
    chromins = []
    rejoined_data = []

    # if there's no data for the primer, then skip to returning empty datasets
    if len(clusterdata) == 0:
        logger.warning(f"No data found for primer '{primer}'.")
        updated_data = clusterdata.copy()
        updated_data["junction"] = ""

    else:
        # drop the seq column for size
        smaller_df = clusterdata.drop(columns=["seq"])

        # group the data by read (qname)
        by_qname = smaller_df.groupby("qname")

        # for each read
        for qname, group in by_qname:
            # check all alignment records for the read have the same value for number of alignments, and
            # get the number of alignments for the read in the group (cluster)
            if len(set(group["n_alignments"])) == 1:
                aln_num = group["n_alignments"].iloc[0]
            else:
                logger.error(f"more than one alignment count for {qname}")
                exit()

            # if user wants to remove inferred alignments (must do this prior to junction calculation)
            if filter_inferred:
                # find all the inferred alignments for the group
                inferred = group[group['inferred_by_primer'] == 1].index
                if not inferred.empty:
                    # recaculate and reset the number of alignments, after the inferred ones are removed
                    aln_num = int(group["n_alignments"].iloc[0]) - len(inferred)
                    group = group.drop(inferred)
                    group["n_alignments"] = aln_num

            # calculate the junction sizes, as a new column
            group['junction'] = group['qend'] - group['qstart'].shift(-1)
            rejoined_data.append(group)

            # add the alignment number to the list (either original, or edited following inferred filter)
            aln_lens.append(aln_num)

            # in fslr output, each group is sorted by q positions
            # get the index of the min and max alignment positions along the read, for each read
            mindex = group["qstart"].idxmin()
            maxdex = group["qstart"].idxmax()
            # if the min and max qstarts are the same*, then set the max index to be the max qend instead
            # *(sometimes alignments overlap, and if there are only two alignments, this results in equal min and max)
            if mindex == maxdex:
                maxdex = group["qend"].idxmax()
            # if the min and max qstarts are the same, and max qend has the same index as the min qstart, log a warning
            if mindex == maxdex and aln_num > 1:
                logger.warning(f"Indexes for min qstart and max qend are the same for {qname}, which has {aln_num} alns")

            # append the min to the list of flanking sizes
            flanking_aln_sizes.append(int(group['aln_size'].loc[mindex]))

            # extract the flanking chromosomes and append to the relevant list
            chromflank.append(group['chrom'].loc[mindex])

            # as long as the min and max are different (i.e. there is more than one alignment), also append the max
            if mindex != maxdex:
                flanking_aln_sizes.append(int(group['aln_size'].loc[maxdex]))
                chromflank.append(group['chrom'].loc[maxdex])

            # if there are more than just the flanking alignments, then fetch those, too
            if aln_num > 2:
                # fetch all alignments from the first non-flanking, to the final non-flanking
                others = group.loc[(mindex+1):(maxdex-1)]  #  in loc, start and the stop are included

                # get aln sizes
                for x in others["aln_size"]:
                    x = int(x)
                    insertion_aln_sizes.append(x)
                    # categorise into simple and complex
                    if aln_num < complex_limit:
                        simple_insert.append(x)
                    else:
                        complex_insert.append(x)

                # get aln chrs
                for x in others["chrom"]:
                    chromins.append(x)

        # make it back into a df again, with the junctions added
        updated_data = pd.concat(rejoined_data)

    return updated_data, aln_lens, flanking_aln_sizes, insertion_aln_sizes, simple_insert, complex_insert, chromflank, chromins

def process_data(clusterdata, primers, complex_limit, logger, keep_false, filter_inferred):
    """
    split the data into required categories
    :param DataFrame clusterdata: df of clustered data
    :param int complex_limit: threshold for defining complexity
    :param list primers: all the primers to include
    :param Logger logger: logger object
    :param bool keep_false: whether to keep alignments with 'false' in the primer suffix
    :param bool filter_inferred: whether to filter inferred alignments
    :return:
    """
    # filter data to specified primers
    filtered = filter_by_primers(clusterdata, primers, logger, keep_false)

    # capture all the required information
    primerfilterdf, aln_lens, flanking_aln_sizes, insertion_aln_sizes, simple_insert, complex_insert, chromflank, chromins = (
        split_flanks_and_inserts(clusterdata=filtered, complex_limit=complex_limit, logger=logger, primer=primers, filter_inferred=filter_inferred))

    # convert lists to dfs for printing to excel
    ins_df = pd.DataFrame({'insertion_aln_sizes':insertion_aln_sizes})
    fla_df = pd.DataFrame({'flanking_aln_sizes':flanking_aln_sizes})
    alns_df = pd.DataFrame({'alns_per_group':aln_lens}).sort_values('alns_per_group', ascending=False)
    simps_df = pd.DataFrame({'simple_insertion_aln_sizes':simple_insert})
    comps_df = pd.DataFrame({'complex_insertion_aln_sizes':complex_insert})

    # count the frequency of each number of alignments
    n_alns_freqs = pd.DataFrame.from_dict(Counter(aln_lens), orient='index').reset_index()
    if n_alns_freqs.empty:  # can be empty, e.g. if primer not present in dataset
        simple_vs_complex = pd.DataFrame()  # create empty df for empty results sheet
    else:
        # rename columns
        n_alns_freqs.columns=['n_alignments', 'count']
        # sort
        n_alns_freqs.sort_values(by='n_alignments', ascending=True, inplace=True)
        # calculate percentages
        n_alns_freqs['percentage'] = n_alns_freqs['count'] / n_alns_freqs['count'].sum() * 100

        # get the combined count of all complex alignments
        all_complex_count = n_alns_freqs[n_alns_freqs['n_alignments']>=complex_limit]['count'].sum()
        # create the df of individual simple alignment counts and total complex alignment counts
        simple = n_alns_freqs[n_alns_freqs['n_alignments']<complex_limit]
        # convert to string, to allow for addition of the string describing the combined complex alns
        simple['n_alignments'] = simple['n_alignments'].astype(str)
        simple_vs_complex = pd.concat([simple, pd.DataFrame({'n_alignments': [f">={str(complex_limit)}"], 'count': [all_complex_count]})], ignore_index=True)
        # calculate percentage
        simple_vs_complex['percentage'] = simple_vs_complex['count'] / simple_vs_complex['count'].sum() * 100

    # count the flanks and inserts, and sort by chromosome (natural sort, for most intuitive order)
    chromflank_counts = pd.DataFrame.from_dict(Counter(chromflank), orient='index').reset_index()
    if not chromflank_counts.empty:
        chromflank_counts.columns=['chrom', 'count']
        chromflank_counts = chromflank_counts.sort_values(by='chrom', key=natsort_keygen())
        # percentage
        chromflank_counts['percentage'] = chromflank_counts['count'] / chromflank_counts['count'].sum() * 100
    chromins_counts = pd.DataFrame.from_dict(Counter(chromins), orient='index').reset_index()
    if not chromins_counts.empty:
        chromins_counts.columns=['chrom', 'count']
        chromins_counts = chromins_counts.sort_values(by='chrom', key=natsort_keygen())
        chromins_counts['percentage'] = chromins_counts['count']/chromins_counts['count'].sum()*100

    # categorise the junctions
    indels = primerfilterdf['junction'][primerfilterdf['junction'] < 0]
    blunts = primerfilterdf['junction'][primerfilterdf['junction'] == 0]
    overhangs = primerfilterdf['junction'][primerfilterdf['junction'] > 0]

    # create a dict with junction counts
    junctions = {'category': ['indels', 'blunts', 'overhangs'],
                 'count': [len(indels), len(blunts), len(overhangs)]}

    # convert to df
    juncs_df = pd.DataFrame(junctions)

    # calculate percentages
    juncs_df['percentage'] = juncs_df['count']/juncs_df['count'].sum()*100

    return primerfilterdf, alns_df, fla_df, ins_df, simps_df, comps_df, n_alns_freqs, chromflank_counts, chromins_counts, simple_vs_complex, juncs_df

def default_output(primerfilterdf, logger, outdir, outfile, type_string, extension):
    """
    if unable to write excel file (e.g. too many rows), write a csv file instead
    :param primerfilterdf: filtered data to be written out
    :param logger: logger object
    :param outdir: path to output directory
    :param outfile: output file name
    :param type_string: type of output file
    :param extension: extension of output file
    :return: tsv file of primerfilterdf
    """
    # create the directory
    res_dir = os.path.join(outdir, "oversized_results")
    os.makedirs(res_dir, exist_ok=True)

    # create and write the output file
    res_tsv = os.path.join(res_dir, f"{os.path.basename(outfile).split(".")[0]}_{type_string}.{extension}")
    primerfilterdf.to_csv(res_tsv, index=False, sep="\t")

    logger.info(f"Done: results written to {res_tsv}")

def run_process(df, primer_strings, complex_limit, logger, orig_name, outdir, outfile, keep_false,output_original, filter_inferred):
    """
    run everything
    :param DataFrame df: dataframe of fslr data to be processed
    :param list primer_strings: list of primer strings to include
    :param int complex_limit: number of alignments to consider as complex
    :param logger logger: logger object
    :param str orig_name: name for the sheet containing the original data
    :param path outdir: output directory
    :param path outfile: path for output file
    :param bool keep_false: whether to keep alignments with false in the primer suffix
    :param bool output_original: whether to output the original data alongside the rest
    :param bool filter_inferred: whether to filter inferred alignments
    :return: excel output file(s)
    """
    # process the data
    primerfilterdf, alns_df, fla_df, ins_df, simps_df, comps_df, n_alns_freqs, chromflank, chromins, simp_vs_comp, junctions = (
        process_data(df, primer_strings, complex_limit, logger, keep_false, filter_inferred))

    # make the full path to the output excel file
    outpath = os.path.join(outdir, outfile)

    default_message = ["Too many rows."]
    default_df = pd.DataFrame(default_message)

    extracted_data = {"simple_insertion_aln_lengths": simps_df,
                      "complex_insertion_aln_lengths": comps_df,
                      "flank_insertion_aln_lengths": fla_df,
                      "alns_per_group": alns_df}

    # write the processed data to the output file, in individual sheets
    logger.info(f"Writing output file '{outpath}'")
    with pd.ExcelWriter(outpath) as writer:
        if output_original:
            logger.info("Writing original data...")
            if len(df) < 1048576:
                # try to write the original data to the first sheet; this may be too large
                try:
                    df.to_excel(writer, sheet_name=orig_name, index=False)
                    logger.info("Done: original data written to excel")
                except ValueError as err:  # if original data is too large, write an almost-empty sheet with an error
                    logger.warning(f"{err}. Skipping to filtered data")
                    default_df.to_excel(writer, sheet_name=orig_name, index=False)
            else:
                logger.warning("Original data too large for Excel. Skipping to filtered data")
                # write an almost-empty first sheet, noting an error
                default_df.to_excel(writer, sheet_name=orig_name, index=False)

        logger.info("Writing filtered data...")
        # if data is expected to be sufficiently small to be written to excel
        if 0 < len(primerfilterdf) < 1048576:
            # try to write the filtered data to additional sheets (should succeed, due to if statement above)
            try:
                primerfilterdf.to_excel(writer, sheet_name="primer_filter", index=False)
                ins_df.to_excel(writer, sheet_name="all_insertion_aln_sizes", index=False)  # same length as filtered
                logger.info("Done: filtered data written to excel")
            except ValueError as err:  # log errors, and output a tsv instead
                logger.warning(f"{err}. Skipping to summary data")
                # write something to the excel file, to prevent IndexError for With Open
                default_df.to_excel(writer, sheet_name="error_encountered", index=False)
        else:
            # if data too large, write to bed instead
            if len(primerfilterdf) > 0:
                logger.warning(f"Output is too large for excel. Producing bed file of filtered results instead.")
                default_output(primerfilterdf, logger, outdir, outfile, "filtered_data", "bed")
                default_output(ins_df, logger, outdir, outfile, "insertion_aln_lengths", "csv")

        # for the variable length sheets, try to write them, log errors, and write to csv if too large for excel
        for k, v in extracted_data.items():
            if 0 < len(v) < 1048576:
                try:
                    v.to_excel(writer, sheet_name=k, index=False)
                except ValueError as err:
                    logger.warning(f"{err}")
            else:  # if data too large, write to csv instead
                if len(v) > 0:
                    default_output(v, logger, outdir, outfile, k, "csv")

        # always at least write out the summary tables to excel
        logger.info("Writing summary data to excel...")
        n_alns_freqs.to_excel(writer, sheet_name="n_alns_freqs", index=False)
        simp_vs_comp.to_excel(writer, sheet_name="combined_complex_freq", index=False)
        chromflank.to_excel(writer, sheet_name="chromflank_counts", index=False)
        chromins.to_excel(writer, sheet_name="chromins_counts", index=False)
        junctions.to_excel(writer, sheet_name="junctions", index=False)
        logger.info("Done: summary data written to excel")

    ### deprecated: writing summaries to excel now, instead
    # # create the summary data output
    # summary_data = (f"Frequency of different numbers of alignments for {primer_strings}:\n"
    #       f"{n_alns_freqs}\n\n"
    #       f"Frequency of simple and complex alignments for {primer_strings}:\n"
    #       f"{simp_vs_comp}\n\n"
    #       f"Frequency of flanking chromosomes for {primer_strings}:\n"
    #       f"{chromflank}\n\n"
    #       f"Frequency of insertion chromosomes for {primer_strings}:\n"
    #       f"{chromins}\n\n"
    #       f"Frequency of different junctions for {primer_strings}:\n"
    #       f"{junctions}")
    #
    # # write the summary to a txt file, even if empty
    # txt_outfile = os.path.join(outdir, f"{os.path.basename(outfile).split(".")[0]}_summary.txt")
    # with open(txt_outfile, "w") as f:
    #     f.write(summary_data)
    #
    # if len(primerfilterdf) > 0:
    #     # log the summary, as long as doing so won't just clog it up with empty tables
    #     logger.info(summary_data)


# CLI setup
@click.command() # Initializes Click for command-line interaction
@click.argument('infiles', nargs=-1, type=click.Path(exists=True))
@click.option('-m', '--mixed_input', is_flag=True, default=False,
              help="Input files should be analysed together (default = analyse each file individually)"
              )
@click.option('-s', '--split', is_flag=True, default=False,
              help="Split output files by primer (i.e. one excel file per provided primer)"
              )
@click.option('-l', '--log_level', type=click.Choice([10, 20, 30, 40, 50]), default=30,
              help='Numeric representation of console logging level to use:'
                   '\n\t10 = DEBUG, 20 = INFO, 30 = WARNING, 40 = ERROR, 50 = CRITICAL (default 30)'
              )
@click.option('-f', '--log_file', is_flag=True, default=False,
              help='Turn on logging to file'
              )
@click.option('-o', '--outdir', type=click.Path(exists=True), default=None,
              help='Path to be used as the output directory (default = same as first input file directory)'
              )
@click.option('-p', '--primer_strings', multiple=True, type=str, default=["21q1"],
              help='primer, as displayed in read names, to use in filtering (default = 21q1)'
              )
@click.option('--prefix', type=click.STRING, default='',
              help='Prefix for output files (default blank)'
              )
@click.option('-c', '--complex_limit', type=int, default=5,
              help='Minimum number of alignments considered to be a complex insertion (default = 5)'
              )
@click.option('-n', '--name', type=click.STRING, default=None,
              help='Base output file name (default = derived from input)'
              )
@click.option('--include_false', is_flag=True, default=False,
              help="Include reads with a false primer (e.g. 'qname1.False_primer1R' or 'qname2.primer1F_False') \n"
                   "(default = False)"
              )
@click.option('--include_original', is_flag=True, default=False,
              help="Print the original data to the first sheet of the output excel file \n"
                   "(only possible if there are fewer rows than the excel file limit, default = False)"
              )
@click.option('--filter_inferred', is_flag=True, default=False,
              help="Filter out alignments inferred by FSLR (and reduce alignment count accordingly, default = False)")

def main(**kwargs):
    """Enter a comma-separated list of input files (FSLR representative.bed files)\n\t"
    "e.g. path/to/my_results1.representative.bed,path/to/my_results2.representative.bed")
    """
    # set the output directory: either as the supplied value, or as the file path of the first input file
    if kwargs['outdir'] is not None:
        outdir = kwargs['outdir']
    else:
        outdir = Path(kwargs['infiles'][0]).parent

    # set up logging, using the output directory determined above
    date_stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    lgr = Logger()
    logger = lgr.initialise_logger(level=kwargs['log_level'], logfile=kwargs['log_file'], outdir=outdir,
                                   prefix=kwargs['prefix'], datestamp=date_stamp)
    logger.info(f"Using {outdir} as the output directory")

    # ensure primer strings unique to prevent duplication of results
    primer_strings = list(set(kwargs['primer_strings']))

    # if the input files need to be combined prior to analysis, then combine them
    if kwargs['mixed_input']:
        # initialise list of all data
        all_dfs = []

        # collate all input data into a single list
        for infile in kwargs['infiles']:
            logger.info(f"Processing {infile}")

            # read file (tsv with header)
            df = pd.read_csv(infile, sep='\t', header=0)
            # record source filename in the output
            df["source"] = os.path.basename(infile).split('.')[0]
            # add the df to the list of dfs
            all_dfs.append(df)

        # concat all the data into a single df
        full_df = pd.concat(all_dfs)

        # set the output file name, either as the supplied value, or as the default
        if kwargs['name'] is not None:
            outfile_name = kwargs['name']
        else:
            # select the file name part of the first input file
            outfile_name = os.path.basename(kwargs['infiles'][0]).split('.')[0]
            if len(kwargs['infiles']) > 1:
                # warn if more than one input file
                logger.warning(f"Creating default output file name from the first of {len(kwargs['infiles'])} files")
        logger.info(f"Using {outfile_name} as the output file name")

        # split the output files, if needed
        if kwargs['split']:
            for ps in primer_strings:
                # generate the full path for the output file
                full_outfile_name = f"{kwargs['prefix']}{ps}_processed_{outfile_name}.xlsx"
                # run the processing, for each primer individually
                run_process(df=full_df, primer_strings=[ps], complex_limit=kwargs['complex_limit'], logger=logger,
                            orig_name='original_data', outdir=outdir, outfile=full_outfile_name,
                            keep_false=kwargs['include_false'], output_original=kwargs['include_original'],
                            filter_inferred=kwargs['filter_inferred'])
        else:
            # generate the full path for the output file
            full_outfile_name = f"{kwargs['prefix']}processed_{outfile_name}.xlsx"
            # run the processing
            run_process(df=full_df, primer_strings=primer_strings, complex_limit=kwargs['complex_limit'], logger=logger,
                        orig_name='original_data', outdir=outdir, outfile=full_outfile_name,
                        keep_false=kwargs['include_false'], output_original=kwargs['include_original'],
                        filter_inferred=kwargs['filter_inferred'])
    else:
        for infile in kwargs['infiles']:
            logger.info(f"Analysing {infile}")
            full_df = pd.read_csv(infile, sep='\t', header=0)

            # set the output file name, either as the supplied value, or as the default
            if kwargs['name'] is not None:
                outfile_name = kwargs['name']
            else:
                # select the file name part of the first input file
                outfile_name = os.path.basename(infile).split('.')[0]

            # split the output files, if needed
            if kwargs['split']:
                for ps in primer_strings:
                    # generate the full path for the output file
                    full_outfile_name = f"{kwargs['prefix']}{ps}_processed_{outfile_name}.xlsx"
                    # run the processing, for each primer individually
                    run_process(df=full_df, primer_strings=[ps], complex_limit=kwargs['complex_limit'], logger=logger,
                                orig_name='original_data', outdir=outdir, outfile=full_outfile_name,
                                keep_false=kwargs['include_false'], output_original=kwargs['include_original'],
                                filter_inferred=kwargs['filter_inferred'])
            else:
                # generate the full path for the output file
                full_outfile_name = f"{kwargs['prefix']}processed_{outfile_name}.xlsx"
                # run the processing
                run_process(df=full_df, primer_strings=primer_strings, complex_limit=kwargs['complex_limit'], logger=logger,
                            orig_name='original_data', outdir=outdir, outfile=full_outfile_name,
                            keep_false=kwargs['include_false'], output_original=kwargs['include_original'],
                            filter_inferred=kwargs['filter_inferred'])

    logger.info(f"Done. Processing complete.")

if __name__ == '__main__':
    main()