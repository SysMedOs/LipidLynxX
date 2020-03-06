# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os
import sys
import pytest

lipidlynx_Path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lipidlynx_Path + "/../")

from lynx.models.log import logger
from lynx.models.rules import InputRules, OutputRules
from lynx.controllers.general_functions import get_abs_path, js_reader

test_input_files = [
    r"../lynx/configurations/rules/input/LipidLynxX.json",
    r"../lynx/configurations/rules/input/MS-DIAL.json",
    r"../lynx/configurations/rules/input/LIPIDMAPS_LMSD.json",
]

test_output_files = [r"../lynx/configurations/rules/output/LipidLynxX.json"]


@pytest.mark.parametrize("test_file", test_input_files)
def test_input_rule(test_file):
    logger.debug("SETUP TESTS...")
    logger.info(test_file)
    # in_file = get_abs_path("test_input/Input_LIPIDMAPS_ShortHand.csv")
    in_file = None
    if test_file:
        in_file = get_abs_path(test_file)
    if not in_file:
        in_file = get_abs_path(r"../lynx/configurations/rules/input/LipidLynxX.json")
    logger.info(f"Test file {in_file}")
    rule = InputRules(test_file)
    logger.debug(f"Got infile {in_file}")
    logger.debug(f"test input rule {rule.source}")
    if rule.is_validated is False:
        raise Exception(f"FAILED: test input rule {rule.source}")
    else:
        logger.info(f"PASSED: test input rule  {rule.source}")
    logger.info(f"test PASSED")


@pytest.mark.parametrize("test_file", test_output_files)
def test_output_rule(test_file):
    logger.debug("SETUP TESTS...")
    logger.info(test_file)
    # in_file = get_abs_path("test_input/Input_LIPIDMAPS_ShortHand.csv")
    in_file = None
    if test_file:
        in_file = get_abs_path(test_file)
    if not in_file:
        in_file = get_abs_path(r"../lynx/configurations/rules/output/LipidLynxX.json")
    logger.info(f"Test file {in_file}")
    rule = OutputRules(test_file)
    logger.debug(f"Got Output infile {in_file}")
    logger.debug(f"test Output rule {rule.nomenclature}")
    if rule.is_structure_valid is False:
        raise Exception(f"FAILED: test Rule {rule.nomenclature}")
    else:
        logger.info(f"PASSED: test Rule {rule.nomenclature}")
        logger.info(f"Supported LMSD classes: {rule.supported_lmsd_classes}")
    logger.info(f"test PASSED")
