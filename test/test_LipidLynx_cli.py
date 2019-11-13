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

lipidlynx_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lipidlynx_Path + "/../")

from lynx_cli import LipidLynx_cli
from lipidlynx.models.DefaultParams import logger


class lipidlynxTestCase(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS...")

        in_file_lst = [
            r"../test/TestInput/test_names.txt",
            r"test/TestInput/test_names.txt",
            r"../TestInput/test_names.txt",
            r"TestInput/test_names.txt",
        ]
        in_file = ""
        for f in in_file_lst:
            if os.path.isfile(f):
                in_file = os.path.abspath(f)
                break
        logger.info(f"Input file is: {in_file}")

        bad_in_file = r"test/TestInput/test_names_x.txt"

        out_folder_lst = [
            r"../test/TestOutput/",
            r"test/TestOutput/",
            r"../TestOutput/",
            r"TestOutput/",
        ]
        out_folder = ""
        for p in out_folder_lst:
            if os.path.isdir(p):
                out_folder = os.path.abspath(p)
                break
        out_file = os.path.join(out_folder, "test_names_sdf.sdf")
        logger.info(f"Out put file will be: {out_file}")

        self.pass_params = ["-i", in_file, "-o", out_file]
        self.fail_input_params = ["-i", bad_in_file, "-o", out_file]

    @staticmethod
    def test_lipidlynx_help():
        logger.debug("test help...")
        result = LipidLynx_cli.main(["-h"])
        if result is False:
            logger.debug("test help... PASSED")
        else:
            raise Exception("test help... Failed")

    @staticmethod
    def test_lipidlynx_bad_params():
        logger.debug("test bad parameter...")
        result = LipidLynx_cli.main(["-test"])
        if result is False:
            logger.debug("test bad parameter... PASSED")
        else:
            raise Exception("test bad parameter... Failed")

    def test_lipidlynx_bad_input(self):
        logger.debug("test bad input...")
        result = LipidLynx_cli.main(self.fail_input_params)
        if result is False:
            logger.debug("test bad input... PASSED")
        else:
            raise Exception("test bad input... Failed")

    def test_lipidlynx_good_input(self):
        logger.debug("test sample data...")
        result = LipidLynx_cli.main(self.pass_params)
        if result is True:
            logger.debug("test sample data... PASSED")
        else:
            raise Exception("test sample data... Failed")

    def tearDown(self):
        logger.debug("TEST END!")


if __name__ == "__main__":
    # python LipidLynx_cli.py -i test/TestInput/test_names.txt -o test/TestOutput/test_names_sdf.sdf

    unittest.main()
    logger.info("TESTS FINISHED!")
