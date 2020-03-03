# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
import os
import sys
import pandas as pd
import pytest

lipidlynx_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lipidlynx_Path + "/../")

from lynx.models.log import logger
from lynx.controllers.rules_reader import Rules, rgx_reader
from lynx.controllers.general_functions import get_abs_path


test_files = [
    r"../lynx/configurations/rules/MS-DIAL.json",
    r"../lynx/configurations/rules/LIPIDMAPS_LMSD.json",
]


@pytest.mark.parametrize("test_file", test_files)
def test_rule(test_file):
    logger.debug("SETUP TESTS...")
    logger.info(test_file)
    # in_file = get_abs_path("test_input/Input_LIPIDMAPS_ShortHand.csv")
    in_file = None
    if test_file:
        in_file = get_abs_path(test_file)
    if not in_file:
        in_file = get_abs_path(r"lynx/configurations/rules/MS-DIAL.json")
    logger.info(f"Test file {in_file}")
    js_dct = rgx_reader(test_file)
    rule = Rules(js_dct)
    logger.debug(f"Got infile {in_file}")
    logger.debug(f"test Rule {rule.source}")
    is_valid = rule.validate()
    if is_valid is False:
        raise Exception(f"FAILED: test Rule {rule.source}")
    else:
        logger.info(f"PASSED: test Rule {rule.source}")
    logger.info(f"test PASSED")
