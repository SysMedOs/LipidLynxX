# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os
import sys
import unittest

from jsonschema import Draft7Validator

epiLION_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, epiLION_Path + "/../")

from lipidlynx.models.DefaultParams import lynx_schema
from lipidlynx.models.DefaultParams import logger
from lipidlynx.controllers.GeneralFunctions import get_abs_path


class JsonTestCase(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS...")
        self.schema_test_dct = {}
        for sn in lynx_schema:
            f = get_abs_path(lynx_schema[sn])
            if f.endswith('.schema.json'):
                n = f'{f[:-12]}.json'
                if os.path.isfile(f) and os.path.isfile(n):
                    self.schema_test_dct[f] = n
                else:
                    raise IOError(f'Can not find schema.json file and json example: {f}')

        if not self.schema_test_dct:
            raise IOError(f'Can not find schema.json file and json example in settings: {json.dumps(lynx_schema)}')

    def test_schema(self):
        logger.debug("test schema...")
        for s in self.schema_test_dct:
            with open(s, 'r') as s_obj:
                s_json = json.load(s_obj)
                with open(self.schema_test_dct[s], 'r') as j_obj:
                    j_json = json.load(j_obj)
                    validator = Draft7Validator(s_json)
                    if validator.is_valid(j_json):
                        logger.debug(f"Schema test PASSED:  {s}")
                    else:
                        raise Exception(f"Schema test FAILED:  {s}")
