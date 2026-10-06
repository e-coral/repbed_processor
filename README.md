# Representative bed file processor
This tool is designed to parse [FLSR](https://github.com/kcleal/fslr) output files (representative bed files), and create excel files containing separate sheets for various categories of data.

## Installation
This is a python tool, available on pypi. To install, run

`pip install repbed_processor`

## Running the tool

### Simple examples

To see more information about the available options for running the tool, run:

`python -m repbed_processor --help`

For a single input file, and using the default probes of interest and default threshold for complex events, simply run:

`python -m repbed_processor path/to/your_file.representative.bed` 

(where 'path/to/your_file.representative.bed' is the path to your FSLR output file)

### Other examples

To change the threshold for what you consider to be a complex event, use `-c` (default is 5):

`python -m repbed_processor path/to/your_file.representative.bed -c 10` 

If you have multiple input files, you may enter them all:

`python -m repbed_processor path/to/your_file1.representative.bed path/to/your_file2.representative.bed`

And if you would like to analyse multiple input files all together (e.g. combining several controls), 
use the `--mixed-input` (or `-m`) flag:

`python -m repbed_processor path/to/your_file1.representative.bed path/to/your_file2.representative.bed --mixed-input`

To specify your probes (rather than using the default, 21q1), use the -p option (you must specify -p before each probe):

`python -m repbed_processor path/to/your_file.representative.bed -p probe_name -p another_probe`

If you would like separate output files for each probe, use the `--split` (or `-s`) flag:

`python -m repbed_processor path/to/your_file.representative.bed -p probe_name -p another_probe -s`

Other options, such as logging and file name customisation, are detailed within the help message (see `--help`, above)

## Issues
Please log issues on [GitHub](https://github.com/e-coral/repbed_processor)

## Licence
MIT (see LICENCE.txt)
