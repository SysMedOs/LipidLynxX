# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from typing import List, Dict, Union, Tuple


from .encoder import lynx_encode
from .parser import parse


def convert_string(
    input_abbr: str, output_dct: Dict[str, Union[List]] = None
) -> Dict[str, Union[List, List[Tuple]]]:
    if output_dct:
        pass
    else:
        output_dct = {"input": [], "output": [], "converted": [], "skipped": []}
    if input_abbr and isinstance(input_abbr, str) and len(input_abbr) < 512:
        lynx_id = lynx_encode(parse(input_abbr))
        if lynx_id:
            output_dct["input"].append(input_abbr)
            output_dct["output"].append(lynx_id)
            output_dct["converted"].append((input_abbr, lynx_id))
        else:
            output_dct["skipped"].append(input_abbr)
    return output_dct


def convert_list(input_lst: List[str]) -> Dict[str, Union[List, List[Tuple]]]:
    output_dct = {"input": [], "output": [], "converted": [], "skipped": []}
    if input_lst and isinstance(input_lst, list):
        input_lst = list(filter(lambda x: isinstance(x, str) and len(x) > 1, input_lst))
        for abbr in input_lst:
            output_dct = convert_string(abbr, output_dct)
    return output_dct


def convert_dict(input_dct: dict) -> Dict[str, Union[List, List[Tuple]]]:
    output_dct = {"input": [], "output": [], "converted": [], "skipped": []}
    if input_dct and isinstance(input_dct, dict):
        for k in input_dct:
            if isinstance(k, str) and len(k) < 256:
                k_val = input_dct[k]
                if isinstance(k_val, list):
                    in_lst = k_val
                else:
                    in_lst = [k_val]
                if in_lst:
                    output_dct[k] = convert_list(in_lst)
    return output_dct
