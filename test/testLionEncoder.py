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

epiLION_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, epiLION_Path + "/../")

from epilion.controllers.Logger import logger
from epilion.controllers.Parser import parse
from epilion.controllers.Encoder import lion_encode
from epilion.controllers.GeneralFunctions import get_abs_path


class ConvertTestCase(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS...")

        in_file = get_abs_path("TestInput/Lipid_Abbreviations.csv")
        self.in_df = pd.read_csv(in_file, header=0, index_col=False)
        logger.debug(f"Got infile {in_file}")

    def test_lion_encode(self):
        logger.debug("test epiLION Encoder ...")
        in_df_test = self.in_df[self.in_df["CONVERT"] == "T"]
        for i, r in in_df_test.iterrows():
            logger.info(f'Process Lipid: {r["INPUT"]}')
            parsed_dct = parse(r["INPUT"])
            test_output = lion_encode(parsed_dct)
            correct_output = r["OUTPUT"].strip('"')
            correct_output = correct_output.strip('"')
            if test_output != correct_output:
                raise Exception(
                    f'FAILED: {r["INPUT"]} -> {test_output} != {correct_output}'
                )
            else:
                logger.info(
                    f'PASSED: {r["INPUT"]} -> {test_output} == {correct_output}'
                )
        logger.info(f"test PASSED")

    @pytest.mark.skip(reason="Not yet finished")
    def test_batch_encode(self):
        logger.debug("test parse_lion ...")
        in_df_test = self.in_df[self.in_df["CONVERT"] == "T"]
        for i, r in in_df_test.iterrows():
            logger.info(f'Process Lipid: {r["INPUT"]}')
            parsed_dct = parse(r["INPUT"])
            test_output = lion_encode(parsed_dct)
            correct_output = r["OUTPUT"].strip('"')
            correct_output = correct_output.strip('"')
