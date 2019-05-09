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

from epilion.controllers import epiLION
from epilion.LibLION.DefaultParams import logger


class epiLIONTestCase(unittest.TestCase):

    def setUp(self):
        logger.debug('SETUP TESTS...')
        in_file = r'test/TestInput/test_names.txt'
        bad_in_file = r'test/TestInput/test_names_x.txt'
        out_file = r'test/TestOutput/test_sdf.sdf'
        self.pass_params = ['-i', in_file, '-o', out_file]
        self.fail_input_params = ['-i', bad_in_file, '-o', out_file]

    @staticmethod
    def test_epiLION_help():
        logger.debug('test help...')
        result = epiLION.main(['-h'])
        if result is False:
            logger.debug('test help... PASSED')
        else:
            raise Exception('test help... Failed')

    @staticmethod
    def test_epiLION_bad_params():
        logger.debug('test bad parameter...')
        result = epiLION.main(['-test'])
        if result is False:
            logger.debug('test bad parameter... PASSED')
        else:
            raise Exception('test bad parameter... Failed')

    def test_epiLION_bad_input(self):
        logger.debug('test bad input...')
        result = epiLION.main(self.fail_input_params)
        if result is False:
            logger.debug('test bad input... PASSED')
        else:
            raise Exception('test bad input... Failed')

    def test_epiLION_good_input(self):
        logger.debug('test sample data...')
        result = epiLION.main(self.pass_params)
        if result is True:
            logger.debug('test sample data... PASSED')
        else:
            raise Exception('test sample data... Failed')

    def tearDown(self):
        logger.debug('TEST END!')


if __name__ == '__main__':
    # python epiLION.py -i test/TestInput/test_names.txt -o test/TestOutput/test_sdf.sdf
    unittest.main()
    logger.info('TESTS FINISHED!')
