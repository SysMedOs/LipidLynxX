# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from typing import List, Dict, Union, Tuple

from lynx.controllers.encoder import Encoder
from lynx.utils.toolbox import keep_string_only


class Converter:
    def __init__(self):
        self.encoder = Encoder()

    def convert_string(
        self, input_str: str, output_dct: Dict[str, Union[List]] = None
    ) -> Dict[str, Union[List, List[Tuple]]]:
        if output_dct:
            pass
        else:
            output_dct = {"input": [], "output": [], "converted": [], "skipped": []}
        if input_str and isinstance(input_str, str) and len(input_str) < 512:
            converted_id = self.encoder.convert(input_str)
            if converted_id:
                output_dct["input"].append(input_str)
                output_dct["output"].append(converted_id)
                output_dct["converted"].append((input_str, converted_id))
            else:
                output_dct["skipped"].append(input_str)
        return output_dct

    def convert_list(
        self, input_list: List[str]
    ) -> Dict[str, Union[List, List[Tuple]]]:
        output_dct = {"input": [], "output": [], "converted": [], "skipped": []}
        if input_list and isinstance(input_list, list):
            input_list = keep_string_only(input_list)
            for abbr in input_list:
                output_dct = self.convert_string(abbr, output_dct)
        return output_dct

    def convert_dict(
        self, input_dct: Dict[str, Union[str, List[str]]]
    ) -> Dict[str, dict]:
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
                        output_dct[k] = self.convert_list(in_lst)

        return output_dct

    def convert(self, data: Union[dict, List[str], str]) -> Dict[str, dict]:
        output_dct = {}

        if isinstance(data, str):
            output_dct = self.convert_string(data)
        elif isinstance(data, list):
            output_dct = self.convert_list(data)
        elif isinstance(data, dict):
            output_dct = self.convert_dict(data)
        else:
            raise TypeError(
                f"Type: {type(data)} not supported. Supported types: str, List[str], and Dict[str, List[str]]."
            )
        return output_dct


if __name__ == "__main__":
    from lynx.utils.log import logger

    t_in_lst = [
        # "GM3(d18:1/18:2(9Z,11Z)(12OH))",
        # "TG P-18:1_18:2(9Z,11Z)(12OH)_18:1(9)(11OH)",
        # "CL(1'-[18:1(9Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",
        "TG(16:0/18:2/9:0<oxo{9}>)",
    ]
    lynx_converter = Converter()
    for t_in in t_in_lst:
        t1_out = lynx_converter.convert(t_in)
        logger.info(f"Input: {t_in} -> Best Output: {t1_out}")

    t2_out = lynx_converter.convert(t_in_lst)
    logger.info(f"Input: {t_in_lst} -> Best Output: {t2_out}")

    t3_out = lynx_converter.convert({"1": t_in_lst})
    logger.info(f"Input: {t_in_lst} -> Best Output: {t3_out}")

    logger.info("fin")
