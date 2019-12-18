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

from jsonschema import Draft7Validator, RefResolver

epiLION_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, epiLION_Path + "/../")

from lynx.models.defaults import lynx_schema_cfg, core_schema, core_schema_path
from lynx.models.defaults import logger
from lynx.controllers.general_functions import get_abs_path


JSON_SCHEMA_INFO = {
    "lynx_mod": {
        "mod_0": r"test/test_input/test_json/mod/lynx_mod_0.json",
        "mod_1": r"test/test_input/test_json/mod/lynx_mod_1.json",
        "mod_2": r"test/test_input/test_json/mod/lynx_mod_2.json",
        "mod_3": r"test/test_input/test_json/mod/lynx_mod_3.json",
        "mod_3.1": r"test/test_input/test_json/mod/lynx_mod_3.1.json",
        "mod_3.2": r"test/test_input/test_json/mod/lynx_mod_3.2.json",
        "mod_4": r"test/test_input/test_json/mod/lynx_mod_4.json",
        "mod_4.1": r"test/test_input/test_json/mod/lynx_mod_4.1.json",
        "mod_4.2": r"test/test_input/test_json/mod/lynx_mod_4.2.json",
    },
    "lynx_fa": {
        "fa_lite_S0": r"test/test_input/test_json/fa/lynx_fa_lite_S0.json",
        "fa_lite_S1": r"test/test_input/test_json/fa/lynx_fa_lite_S1.json",
        "fa_lite_S2": r"test/test_input/test_json/fa/lynx_fa_lite_S2.json",
        "fa_lite_S3": r"test/test_input/test_json/fa/lynx_fa_lite_S3.json",
        "fa_lite_S3.1": r"test/test_input/test_json/fa/lynx_fa_lite_S3.1.json",
        "fa_lite_S3.2": r"test/test_input/test_json/fa/lynx_fa_lite_S3.2.json",
        "fa_lite_S4": r"test/test_input/test_json/fa/lynx_fa_lite_S4.json",
        "fa_lite_S4.1": r"test/test_input/test_json/fa/lynx_fa_lite_S4.1.json",
        "fa_lite_S4.2": r"test/test_input/test_json/fa/lynx_fa_lite_S4.2.json",
    },
    # "lynx_core": {
    #     "core_S4.1": r"test/test_input/test_json/lynx_core_S4.1.json",
    #     "core_S4.2": r"test/test_input/test_json/lynx_core_S4.2.json",
    # },
}


class JsonTestCase(unittest.TestCase):
    def setUp(self):
        logger.debug("SETUP TESTS...")
        self.schema_test_dct = {}
        for s_name in lynx_schema_cfg:
            if s_name != "_version":
                f = get_abs_path(lynx_schema_cfg[s_name])
                if f.endswith(".schema.json"):
                    if os.path.isfile(f):
                        self.schema_test_dct[s_name] = f
                    else:
                        raise IOError(
                            f"Cannot find schema.json file and json example: {f}"
                        )
            else:
                logger.info(f"API schema version: {lynx_schema_cfg[s_name]}")

        if not self.schema_test_dct:
            raise IOError(
                f"Cannot find schema.json file and json example in settings: {json.dumps(lynx_schema_cfg)}"
            )

    def test_schema(self):
        logger.debug("Start to test schema...")
        for s_name in self.schema_test_dct:
            json_file_dct = JSON_SCHEMA_INFO[s_name]
            for s_v in json_file_dct:
                j_name = json_file_dct[s_v]
                j_f_path = get_abs_path(j_name)
                with open(self.schema_test_dct[s_name], "r") as s_obj:
                    s_json = json.load(s_obj)
                    with open(j_f_path, "r") as j_obj:
                        j_json = json.load(j_obj)
                        validator = Draft7Validator(
                            s_json,
                            resolver=RefResolver(
                                f"file://{core_schema_path}", referrer=core_schema
                            ),
                        )
                        logger.info(f"Constructed... file://{core_schema_path}")
                        if validator.is_valid(j_json):
                            logger.info(f"Schema test PASSED:  {s_v}")
                            logger.debug(f"json file tested: {j_f_path}")
                        else:
                            errors = validator.iter_errors(j_json)
                            for e in errors:
                                logger.error(e)
                            raise Exception(
                                f"Schema test FAILED: {s_v} \n"
                                f"JSON file: {j_f_path}\n"
                                f"Schema: {s_name}\n"
                                f"Schema file: {self.schema_test_dct[s_name]}"
                            )
