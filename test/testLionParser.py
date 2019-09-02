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

epiLION_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, epiLION_Path + "/../")

from epilion.controllers.Logger import logger
from epilion.controllers.Parser import parse
from epilion.controllers.LionEncoder import lion_encode
from epilion.controllers.GeneralFunctions import get_abs_path


class ConvertTestCase(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS...")

        in_file = get_abs_path("TestInput/FA.csv")
        self.in_df = pd.read_csv(in_file, names=["INPUT", "OUTPUT"], index_col=False)
        logger.debug(f"Got infile {in_file}")

    # def test_parse_epilion(self):
    #     logger.debug("test parse_epilion...")
    #     for i, r in self.in_df.iterrows():
    #         test_output = parse_epilion(r["INPUT"])["id"]
    #         correct_output = r["OUTPUT"].strip('"')
    #         correct_output = correct_output.strip('"')
    #         if test_output != correct_output:
    #             raise Exception(
    #                 f'input: {r["INPUT"]} -> {test_output} != output: {correct_output}'
    #             )
    #         else:
    #             logger.info(
    #                 f'input: {r["INPUT"]} -> {test_output} == output: {correct_output}'
    #             )

    def test_lion_encode(self):
        logger.debug("test parse_lion ...")
        for i, r in self.in_df.iterrows():
            parsed_dct = parse(r["INPUT"])
            test_output = lion_encode(parsed_dct)
            correct_output = r["OUTPUT"].strip('"')
            correct_output = correct_output.strip('"')
            if test_output != correct_output:
                raise Exception(
                    f'input: {r["INPUT"]} -> {test_output} '
                    f"!= output: {correct_output}"
                )
            else:
                logger.info(
                    f'input: {r["INPUT"]} -> {test_output} == output: {correct_output}'
                )
