# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
import os
import sys
import unittest
import pandas as pd
import pytest

lipidlynx_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lipidlynx_Path + "/../")

from lipidlynx.controllers.Logger import logger
from lipidlynx.controllers.Parser import parse
from lipidlynx.controllers.Encoder import lynx_encode
from lipidlynx.controllers.GeneralFunctions import get_abs_path


test_files = [
    r"test/TestInput/Input_General.csv",
    r"test/TestInput/Input_ALEX.csv",
    r"test/TestInput/Input_LIPIDMAPS_ShortHand.csv",
    r"test/TestInput/Input_LipidMatch.csv",
    r"test/TestInput/Input_LPPtiger.csv",
]


@pytest.mark.parametrize("test_file", test_files)
def test_lion_encode(test_file):
    logger.debug("SETUP TESTS...")
    logger.info(test_file)
    # in_file = get_abs_path("TestInput/Input_LIPIDMAPS_ShortHand.csv")
    in_file = None
    if test_file:
        in_file = get_abs_path(test_file)
    if not in_file:
        in_file = get_abs_path("TestInput/Input_LIPIDMAPS_ShortHand.csv")
    logger.info(f"Test file {in_file}")
    in_df = pd.read_csv(in_file, header=0, index_col=False)
    logger.debug(f"Got infile {in_file}")

    logger.debug("test epiLION Encoder ...")
    in_df_test = in_df[in_df["CONVERT"] == "T"]
    for i, r in in_df_test.iterrows():
        logger.info(f'Process Lipid: {r["INPUT"]}')
        parsed_dct = parse(r["INPUT"])
        test_output = lynx_encode(parsed_dct)
        correct_output = r["OUTPUT"].strip('"')
        correct_output = correct_output.strip('"')
        if test_output != correct_output:
            raise Exception(
                f'FAILED: {r["INPUT"]} -> {test_output} != {correct_output}'
            )
        else:
            logger.info(f'PASSED: {r["INPUT"]} -> {test_output} == {correct_output}')
    logger.info(f"test PASSED")
