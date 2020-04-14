# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import getopt
import sys

from lynx.controllers.converter import Converter
from lynx.utils.file_handler import create_converter_output, get_abs_path, get_table
from lynx.utils.log import logger
from lynx.utils.toolbox import keep_string_only


def main(argv):
    """
    :param argv: -i <input LipidLynxX abbreviation file in .csv/.xlsx format> -o <output .csv/.xlsx file>
    """

    in_file = ""
    out_file = ""

    is_output = False
    abs_out_file = ""

    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["help", "infile=", "outfile="])
        logger.debug(f"User input: {opts}, {args}")
    except getopt.GetoptError:
        logger.info(
            "-i <input .csv/.xlsx file> -o <output .csv/.xlsx file>"
        )
        return is_output
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            logger.info(
                "-i <input .csv/.xlsx file> -o <output .csv/.xlsx file>"
            )
            return is_output
        elif opt in ("-i", "--infile"):
            in_file = arg
        elif opt in ("-o", "--outfile"):
            out_file = arg

    if in_file and out_file:
        table_info = get_table(in_file)
        logger.info(f"Load input file: {in_file}")
        table_dct = keep_string_only(table_info)
        if table_dct:
            converter = Converter()
            converted_dct = converter.convert_dict(table_dct)
            abs_out_file = create_converter_output(converted_dct, output_name=out_file)
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
