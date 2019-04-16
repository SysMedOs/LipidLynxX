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

epiLION_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, epiLION_Path + '/../')

import convLION
from LibLION.DefaultParams import logger


class epiLION_ConverterTestCase(unittest.TestCase):

    def setUp(self):
        logger.debug('SETUP TESTS...')
        in_file = r'Test/TestInput/test_crosscheck.xlsx'
        bad_in_file = r'Test/TestInput/test_crosscheck_x.txt'
        out_file = r'Test/TestOutput/test_crosscheck_output.xlsx'
        self.pass_params = ['-i', in_file, '-o', out_file]
        self.fail_input_params = ['-i', bad_in_file, '-o', out_file]

    def test_epiLION_Converter_help(self):
        logger.debug('Test help...')
        assert convLION.main(['-h']) is False

    def test_epiLION_Converter_bad_params(self):
        logger.debug('Test bad params...')
        assert convLION.main(['-test']) is False

    def test_epiLION_Converter_bad_input(self):
        logger.debug('Test bad input...')
        assert convLION.main(self.fail_input_params) is False

    def test_epiLION_Converter_good_input(self):
        logger.debug('Test sample data...')
        assert convLION.main(self.pass_params) is True

    def tearDown(self):
        logger.debug('TEST END!')


if __name__ == '__main__':
    # python convLION.py -i Test/TestInput/test_crosscheck.xlsx -o Test/TestOutput/test_crosscheck_output.xlsx
    unittest.main()
    logger.info('TESTS FINISHED!')
