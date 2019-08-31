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

import epiLION
from epilion.libLION.DefaultParams import logger
from epilion.controllers.LionParser import parse_epilion


class ConvertTestCase(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS...")

        in_file_lst = [
            r"../test/TestInput/FA.csv",
            r"test/TestInput/FA.csv",
            r"../TestInput/FA.csv",
            r"TestInput/FA.csv",
        ]
        in_file = ""
        for f in in_file_lst:
            if os.path.isfile(f):
                in_file = os.path.abspath(f)
                break
        logger.debug("Got infile {in_file}")
        self.in_df = pd.read_csv(in_file, names=["INPUT", "OUTPUT"], index_col=False)
        logger.info(f"Input file is: {in_file}")

    def test_parse_epilion(self):
        logger.debug("test help...")
        for i, r in self.in_df.iterrows():
            test_output = parse_epilion(r["INPUT"])['id']
            if test_output != r["OUTPUT"]:
                raise Exception(f'input: {r["INPUT"]} -> {test_output} '
                                f'!= output: {r["OUTPUT"]}')
