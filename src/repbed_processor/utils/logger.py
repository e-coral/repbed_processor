# create and initialise the logger, based on https://github.com/ecoral/LitSpy/blob/master/litspy/logger.py"""

import pathlib
import logging


class Logger:
    """
    class for creating and initialising the logger object
    """
    def __init__(self):
        """initialise the dict of logging levels"""
        self.log_levels_dict = {10: 'DEBUG', 20: 'INFO', 30: 'WARNING', 40: 'ERROR', 50: 'CRITICAL'}

    @staticmethod
    def get_logfile_name(outdir, name, datestamp):
        """
        initialise the log file name
        :param Path outdir: path to output directory
        :param name: string name for the log file
        :param datestamp: string representing date and time the process started
        :return: Path object for log file
        """
        # initialise the additional info to include in the log file name
        filename_additions = []

        # capture the additional filename information in a string
        if name:
            filename_additions.append(name)
        if datestamp:
            filename_additions.append(datestamp)

        fname_addition = "_".join(filename_additions)

        # if there are any filename additions supplied, add them to the filename
        if len(fname_addition) > 1:
            file_path = pathlib.Path(f'{outdir}/log_{fname_addition}.log')
        else:
            file_path = pathlib.Path(f'{outdir}/log.log')
            
        return file_path

    def initialise_logging_to_file(self, logger, outdir, name, datestamp):
        """
        initialise logging to a file with a file name (optional datestamp, filename), and log formatting

        :param logging.Logger logger: logger object
        :param Path outdir: path to output directory
        :param str name: string name for the log file
        :param str datestamp: string representing date and time the process started
        :return: logger object with configured file handler
        :rtype: logging.Logger
        """
        file_path = self.get_logfile_name(outdir, name, datestamp)

        # create the file handler
        file_h = logging.FileHandler(file_path, encoding='utf-8')

        # create and apply the formatter to the file handler
        file_fmt = logging.Formatter('%(asctime)s \t%(levelname)s:\t%(message)s', datefmt='%Y/%m/%d %H:%M:%S')
        file_h.setFormatter(file_fmt)

        # set file logging level to info
        file_h.setLevel(logging.INFO)

        # add the file handler to the logger object
        logger.addHandler(file_h)
        return logger

    @staticmethod
    def initialise_logging_to_console(int_level, logger):
        """
        initialise logging at the console level

        :param int int_level: logging level (int)
        :param logging.Logger logger: logger object
        :return: logger object with configured console handler
        :rtype: logging.Logger
        """
        # initialise console handler and set the console logging level
        cons_h = logging.StreamHandler()
        cons_h.setLevel(int_level)

        # set console formatter and add it to the console handler
        cons_fmt = logging.Formatter('%(asctime)s|%(levelname)s: %(message)s', datefmt='%H:%M:%S')
        cons_h.setFormatter(cons_fmt)

        # add the handler to the logger and return the logger
        logger.addHandler(cons_h)
        return logger

    def initialise_logger(self, level, logfile, outdir, prefix, datestamp):
        """
        initialise the logger to be used

        :param int level: supplied console logging level
        :param bool logfile: whether to log to file
        :param path outdir: path to output directory
        :param str prefix: supplied filenames prefix
        :param str datestamp: string representing date and time the process started
        :return: logger object
        :rtype: logging.Logger
        """
        # initialise the logger
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)

        # if logging to file has been specified, initialise logging to file
        if logfile:
            logger = self.initialise_logging_to_file(logger, outdir, prefix, datestamp)

        # get and set the console logging level
        logger = self.initialise_logging_to_console(level, logger)

        # log relevant messages about logging initialisation
        logger.info(f"Console logging level set to '{self.log_levels_dict[level]}'")
        if not logfile:
            logger.info(f"Logging to console only. To turn on logging to file, run the command again with the -f flag")

        return logger

