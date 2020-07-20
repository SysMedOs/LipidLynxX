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

import re
from typing import List, Dict, Union, Tuple

from lynx.controllers.encoder import Encoder

from lynx.models.api_models import ConvertedStrData, ConvertedListData
from lynx.utils.log import app_logger
from lynx.models.defaults import (
    default_output_rules,
    default_input_rules,
)
from lynx.utils.toolbox import keep_string_only


class Converter:
    def __init__(
        self,
        style: str = "LipidLynxX",
        input_rules: dict = default_input_rules,
        output_rules: dict = default_output_rules,
        logger=app_logger,
    ):
        self.style = style
        self.encoder = Encoder(
            style=style,
            input_rules=input_rules,
            output_rules=output_rules,
            logger=logger,
        )
        self.logger = logger

    def convert_str(self, input_str: str, level: str = None,) -> ConvertedStrData:
        output_dct = {}
        # Set COMP_DB to max level B2
        if re.search(r"COMP\\s*[_]?\\s*(DB)?", self.style):
            level = "B2"
        if input_str and isinstance(input_str, str) and len(input_str) < 512:
            converted_id = self.encoder.convert(input_str, level=level)
            if converted_id:
                output_dct["input"] = input_str
                output_dct["output"] = converted_id
                output_dct["converted"] = (input_str, converted_id)
            else:
                output_dct["skipped"] = input_str
        converted_str_obj = ConvertedStrData(
            input=output_dct.get("input", input_str),
            output=output_dct.get("output", ""),
            converted=output_dct.get("converted", tuple()),
            skipped=output_dct.get("skipped", ""),
        )
        return converted_str_obj

    def convert_list(
        self, input_list: List[str], level: str = None
    ) -> ConvertedListData:
        output_dct = {"input": [], "output": [], "converted": [], "skipped": []}
        abbr_result_lst = []
        if input_list and isinstance(input_list, list):
            input_list = keep_string_only(input_list, self.logger)
            for abbr in input_list:
                abbr_result_lst.append(self.convert_str(abbr, level=level).dict())
        for abbr_result in abbr_result_lst:
            for k in output_dct:
                output_dct[k].append(abbr_result.get(k, ""))
        for k in output_dct:
            output_dct[k] = [v for v in output_dct[k] if v]  # remove "" or None
        converted_lst_obj = ConvertedListData(
            input=output_dct.get("input"),
            output=output_dct.get("output"),
            converted=output_dct.get("converted"),
            skipped=output_dct.get("skipped"),
        )

        return converted_lst_obj

    def convert_dict(
        self, input_dct: Dict[str, Union[str, List[str]]], level: str = None
    ) -> Dict[str, ConvertedListData]:
        output_dct = {}
        if input_dct and isinstance(input_dct, dict):
            for k in input_dct:
                if isinstance(k, str) and len(k) < 256:
                    k_val = input_dct[k]
                    if isinstance(k_val, list):
                        in_lst = k_val
                    else:
                        in_lst = [k_val]
                    if in_lst:
                        output_dct[k] = self.convert_list(in_lst, level=level)

        return output_dct

    def convert(
        self, data: Union[dict, List[str], str], level: str = None
    ) -> Dict[str, dict]:
        output_dct = {}

        if isinstance(data, str):
            output_dct = self.convert_str(data, level=level)
        elif isinstance(data, list):
            output_dct = self.convert_list(data, level=level)
        elif isinstance(data, dict):
            output_dct = self.convert_dict(data, level=level)
        else:
            raise TypeError(
                f"Type: {type(data)} not supported. Supported types: str, List[str], and Dict[str, List[str]]."
            )
        return output_dct


if __name__ == "__main__":
    # from lynx.utils.log import logger

    t_in_lst = [
        # "GM3(d18:1/18:2(9Z,11Z)(12OH))",
        # "TG P-18:1_18:2(9Z,11Z)(12OH)_18:1(9)(11OH)",
        # "TG P-18:1_18:2(9Z,11Z)_18:1(9)",
        # "CL(1'-[18:1(9Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",
        # "TG(16:0/18:2/9:0<oxo{9}>)",
        # "HETE",
        # "HETE",
        # "SPBP 18:0;O",
        # "SPBP 18:0;O3",
        # "Cer 18:1;3O/20:4",
        # "CoA(20:3(11Z,14Z,17Z))",
        # "CoA 18:2;O",
        # "FACoA 18:0",
        # "Cer 24:2",
        "LMGP01010594",
    ]
    lv = "B1"
    # test_out_rule = "COMP_DB"
    test_out_rule = "LipidLynxX"
    lynx_converter = Converter(style=test_out_rule, logger=app_logger)
    for t_in in t_in_lst:
        t1_out = lynx_converter.convert(t_in, level="B1")
        app_logger.info(f"Input: {t_in} -> Best Output: {t1_out}")

    # t2_out = lynx_converter.convert(t_in_lst)
    # logger.info(f"Input: {t_in_lst} -> Best Output: {t2_out}")
    #
    # t3_out = lynx_converter.convert({"1": t_in_lst})
    # logger.info(f"Input: {t_in_lst} -> Best Output: {t3_out}")

    # logger.info("fin")
