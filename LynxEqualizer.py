# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import getopt
import sys

from lynx.controllers.equalizer import Equalizer
from lynx.utils.file_handler import create_equalizer_output, get_abs_path, get_table
from lynx.utils.log import logger
from lynx.utils.toolbox import keep_string_only


def main(argv):
    """
    :param argv: -l<list of levels separated by , e.g.: "B0,D0,D1">
        -i <input LipidLynxX abbreviation file in .csv/.xlsx format> -o <output .csv/.xlsx file>
    """

    levels = []
    in_file = ""
    out_file = ""
    equalized_dct = {}

    is_output = False
    abs_out_file = ""

    try:
        opts, args = getopt.getopt(argv, "hl:i:o:", ["help", "levels=", "infile=", "outfile="])
        logger.debug(f"User input: {opts}, {args}")
    except getopt.GetoptError:
        logger.info(
            "-i <input .csv/.xlsx file> -o <output .csv/.xlsx file>"
        )
        return is_output
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            logger.info(
                '-l<list of levels separated by , e.g.: "B0,D0,D1"> '
                '-i <input .csv/.xlsx file> -o <output .csv/.xlsx file>'
            )
            return is_output
        elif opt in ("-l", "--levels"):
            arg = arg.strip("'")
            arg = arg.strip('"')
            levels = arg.split(',')
            levels = [_lv.strip(" ") for _lv in levels]
        elif opt in ("-i", "--infile"):
            in_file = arg
        elif opt in ("-o", "--outfile"):
            out_file = arg
    if in_file:
        table_info = get_table(in_file)
        logger.info(f"Load input file: {in_file}")
        table_dct = keep_string_only(table_info)
        if table_dct:
            for lv in levels:
                equalizer = Equalizer(input_data=table_dct, level=lv)
                equalized_dct[lv] = equalizer.cross_match()

            abs_out_file = create_equalizer_output(equalized_dct, output_name=out_file)
            if abs_out_file.lower().endswith('.xlsx'):
                abs_out_file = get_abs_path(abs_out_file)
                is_output = True

    if is_output:
        logger.info(f"Save output file: {abs_out_file}")
        logger.info("FINISHED")

        logger.info(f"is_output {is_output}")
    else:
        logger.error(f"Can NOT open input file:")
        logger.error(in_file)
        logger.error("!!! FAILED to PROCESS !!!")

    return is_output


if __name__ == "__main__":

    main(sys.argv[1:])
