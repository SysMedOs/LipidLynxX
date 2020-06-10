# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# LipidLynxX is Dual-licensed
#   For academic and non-commercial use: GPLv2 License:
#   For commercial use: please contact the SysMedOs team by email.
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: lipid annotations converter for large scale lipidomics and epilipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os
import sys
import pytest

lynx_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lynx_path + "/../")

from lynx.utils.log import logger
from lynx.models.rules import InputRules, OutputRules
from lynx.utils.file_handler import get_abs_path

test_input_files = [
    r"lynx/configurations/rules/input/LipidLynxX.json",
    r"/lynx/configurations/rules/input/ALEX123.json",
    r"lynx/configurations/rules/input/LDA.json",
    r"lynx/configurations/rules/input/LipidBlast.json",
    r"lynx/configurations/rules/input/LipiDex.json",
    r"lynx/configurations/rules/input/LipidHome.json",
    r"lynx/configurations/rules/input/LipidLynxX.json",
    r"lynx/configurations/rules/input/LIPIDMAPS_LMSD.json",
    r"lynx/configurations/rules/input/LipidMatch.json",
    r"lynx/configurations/rules/input/LipidPro.json",
    r"lynx/configurations/rules/input/LPPtiger.json",
    r"lynx/configurations/rules/input/MS-DIAL.json",
    r"lynx/configurations/rules/input/Shorthand_based.json",
]

test_output_files = [r"lynx/configurations/rules/output/LipidLynxX.json"]


@pytest.mark.parametrize("test_file", test_input_files)
def test_input_rule(test_file):
    logger.debug("SETUP TESTS...")
    logger.info(test_file)
    in_file = None
    if test_file:
        in_file = get_abs_path(test_file)
    if not in_file:
        in_file = get_abs_path(r"lynx/configurations/rules/input/LipidLynxX.json")
    logger.info(f"Test file {in_file}")
    rule = InputRules(test_file)
    logger.debug(f"Got infile {in_file}")
    logger.debug(f"test input rule {rule.sources}")
    if rule.is_validated is False:
        raise Exception(f"FAILED: test input rule {rule.sources}")
    else:
        logger.info(f"PASSED: test input rule  {rule.sources}")
    logger.info(f"test PASSED")


@pytest.mark.parametrize("test_file", test_output_files)
def test_output_rule(test_file):
    logger.debug("SETUP TESTS...")
    logger.info(test_file)
    in_file = None
    if test_file:
        in_file = get_abs_path(test_file)
    if not in_file:
        in_file = get_abs_path(r"lynx/configurations/rules/output/LipidLynxX.json")
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
