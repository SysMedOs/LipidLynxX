# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import getopt
import os.path
import sys

from LibLION.DefaultParams import logger
from LibLION.Converter import Converter


# required to perform multiprocessing
# import multiprocessing


def main(argv):
    """
    :param argv: -i <input epiLION abbreviation file in .txt format>
    """

    in_file = ''
    out_file = ''

    is_output = False

    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["infile=", "outfile="])
        logger.debug(f'User input: {opts}, {args}')
    except getopt.GetoptError:
        logger.info('epiLIONConverter.py -i <input_file> -o <output_file>')
        return is_output
    for opt, arg in opts:
        if opt == '-h':
            logger.info('epiLIONConverter.py -i <input_file> -o <output_file>')
            return is_output
        elif opt in ("-i", "--infile"):
            in_file = arg
        elif opt in ("-o", "--outfile"):
            out_file = arg

    if os.path.isfile(in_file):
        logger.info(f'Load input file: {in_file}')
        cfg_file = r'Configurations/LinearFA_abbreviations.xlsx'
        if os.path.isfile(cfg_file):
            pass
        elif os.path.isfile(f'../{cfg_file}'):
            cfg_file = f'../{cfg_file}'
        converter = Converter(cfg_file)
        converter.convert_table(in_file, out_file)

        logger.info(f'Save output file: {out_file}')
        logger.info('FINISHED')
        is_output = True
        logger.info(f'is_output {is_output}')
    else:
        logger.error(f'Can NOT open input file:')
        logger.error(in_file)
        logger.error('!!! FAILED to PROCESS !!!')

    return is_output


if __name__ == "__main__":
    # multiprocessing.freeze_support()
    main(sys.argv[1:])
