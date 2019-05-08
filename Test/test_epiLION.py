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

import epiLION
from LibLION.DefaultParams import logger


class epiLIONTestCase(unittest.TestCase):

    def setUp(self):
        logger.debug('SETUP TESTS...')
        in_file = r'Test/TestInput/test_names.txt'
        bad_in_file = r'Test/TestInput/test_names_x.txt'
        out_file = r'Test/TestOutput/test_sdf.sdf'
        self.pass_params = ['-i', in_file, '-o', out_file]
        self.fail_input_params = ['-i', bad_in_file, '-o', out_file]

    @staticmethod
    def test_epiLION_help():
        logger.debug('Test help...')
        result = epiLION.main(['-h'])
        if result is False:
            logger.debug('Test help... PASSED')
        else:
            raise Exception('Test help... Failed')

    @staticmethod
    def test_epiLION_bad_params():
        logger.debug('Test bad parameter...')
        result = epiLION.main(['-test'])
        if result is False:
            logger.debug('Test bad parameter... PASSED')
        else:
            raise Exception('Test bad parameter... Failed')

    def test_epiLION_bad_input(self):
        logger.debug('Test bad input...')
        result = epiLION.main(self.fail_input_params)
        if result is False:
            logger.debug('Test bad input... PASSED')
        else:
            raise Exception('Test bad input... Failed')

    def test_epiLION_good_input(self):
        logger.debug('Test sample data...')
        result = epiLION.main(self.pass_params)
        if result is True:
            logger.debug('Test sample data... PASSED')
        else:
            raise Exception('Test sample data... Failed')

    def tearDown(self):
        logger.debug('TEST END!')


if __name__ == '__main__':
    # python epiLION.py -i Test/TestInput/test_names.txt -o Test/TestOutput/test_sdf.sdf
    unittest.main()
    logger.info('TESTS FINISHED!')
