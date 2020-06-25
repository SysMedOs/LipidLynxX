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

from natsort import natsorted
import pandas as pd
import pytest

lynx_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lynx_path + "/../")

from lynx.controllers.converter import Converter

from lynx.models.defaults import supported_levels
from lynx.utils.cfg_reader import app_cfg_info
from lynx.utils.file_handler import get_abs_path
from lynx.utils.log import app_logger
from lynx.utils.params_loader import (
    build_input_rules,
    build_output_rules,
)


default_input_rules = build_input_rules(app_cfg_info["input_rules"], app_logger)
default_output_rules = build_output_rules(app_cfg_info["output_rules"], app_logger)

default_test_lipids = [
    # ("PLPC", "B1", "LipidLynxX", "PC(34:2)"),
    # ("Cer 24:2", "S1", "LipidLynxX", "Cer(18:1;O2/24:2)"),
    # ("Cer 24:2", "B2", "COMP_DB", "Cer 42:3;O2"),
    # ("dhCer 16:0", "B2", "COMP_DB", "Cer 34:0;O2"),
    # ("DG(O-16:0/18:1)", "B2", "COMP_DB", "DG O-34:1"),
    # ("PC(O-32:1)", "B2", "COMP_DB", "PC O-32:1"),
    # ("CerP 24:2", "B2", "COMP_DB", "CerP 42:3;O2"),
    ("FA 18:2(9,11);O", "B2", "COMP_DB", "FA 18:2;O"),
    # ("CerP 24:2", "S1", "LipidLynxX", "CerP(18:1;O2/24:2)"),
]


default_test_files = [
    ("COMP_DB", r"test/test_input/test_style_compdb.csv"),
    # ("ShorthandNotation", r"test/test_input/test_style_shorthand.csv"),
]


@pytest.mark.parametrize("lipid,level,style,converted_lipid", default_test_lipids)
def test_convert_lipid(
    lipid: str, level: str, style: str, converted_lipid: str,
):
    print(
        f"Convert {lipid} into {level} Level using {style} Style as {converted_lipid}."
    )
    lynx_converter = Converter(
        style=style,
        input_rules=default_input_rules,
        output_rules=default_output_rules,
        logger=app_logger,
    )
    converted_name = lynx_converter.convert_str(input_str=lipid, level=level).output
    assert converted_name == converted_lipid


@pytest.mark.parametrize("style,file", default_test_files)
def test_convert_file(
    style: str, file: str,
):
    in_df = pd.read_csv(get_abs_path(file))
    in_df.fillna("", inplace=True)
    test_df = pd.DataFrame(data=in_df[in_df["SUPPORTED"] != ""])  # type: pd.DataFrame
    levels = test_df.columns.values.tolist()
    levels = [t_lv for t_lv in levels if t_lv in supported_levels]
    max_level = natsorted(levels)[-1]
    test_df.loc[:, "MAX"] = test_df[max_level]
    all_levels = supported_levels.copy()
    all_levels.append("MAX")
    sum_test_dct = test_df.to_dict(orient="index")
    for idx in sum_test_dct:
        test_dct = sum_test_dct[idx]
        for lv in all_levels:
            if lv in test_dct:
                compatible_levels = test_dct.get("SUPPORTED", "T")
                if compatible_levels == "T":
                    lv_compatible_levels = all_levels[
                        : max(all_levels.index(lv) + 1, len(all_levels))
                    ]
                else:
                    lv_compatible_levels = compatible_levels.split(",")
                    lv_compatible_levels = [c_lv.strip() for c_lv in lv_compatible_levels]
                    lv_compatible_levels = [c_lv for c_lv in lv_compatible_levels if c_lv in all_levels]
                test_input = test_dct.get(lv, None)
                for out_lv in lv_compatible_levels:
                    test_output = test_dct.get(out_lv, None)
                    if test_output and test_input:
                        test_convert_lipid(
                            lipid=test_input,
                            level=out_lv,
                            style=style,
                            converted_lipid=test_output,
                        )
