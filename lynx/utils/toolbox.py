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

import base64
import json
import re
from typing import Dict, List, Union

from jsonschema import Draft7Validator

from lynx.models.api_models import level_rgx_str, LvType, StyleType
from lynx.utils.log import logger


def seg_to_str(in_list: List[str], sep: str = ",") -> str:
    """
    Combine multiple segments into one str without space and separator in the end
    Args:
        in_list: input list to be joined.
        sep: separator used to join str segments, default value is ","

    Returns:
        out_str: the joined str

    """
    in_list = filter(None, in_list)
    out_str = sep.join(in_list)
    out_str = out_str.strip(sep)
    if out_str is None:
        out_str = ""
    return out_str


def check_json(validator: Draft7Validator, json_obj: json) -> bool:
    is_valid = False
    if validator.is_valid(json_obj):
        logger.debug(f"JSON Schema check PASSED.")
        is_valid = True
    else:
        for e in validator.iter_errors(json_obj):
            logger.error(e)
        raise Exception(f"JSON Schema check FAILED. Json string: {json_obj}")
    return is_valid


def clean_dct(dct: dict) -> dict:

    if dct:
        for k in dct:
            v = dct[k]
            if isinstance(v, str):
                dct[k] = [v]
            elif isinstance(v, list):
                dct[k] = list(filter(None, v))
            else:
                pass
    else:
        pass

    return dct


def keep_string_only(
    data: Union[list, Dict[str, list]]
) -> Union[list, Dict[str, list]]:
    filtered_data = None
    if isinstance(data, list):
        filtered_data = [d for d in data if isinstance(d, str) and len(d) > 0]
    elif isinstance(data, dict):
        filtered_data = {}
        for k in data:
            filtered_data[k] = [v for v in data[k] if isinstance(v, str) and len(v) > 0]
    else:
        raise TypeError

    if filtered_data:
        return filtered_data
    else:
        raise ValueError


def get_level(lv) -> str:
    if isinstance(lv, LvType):
        use_level = lv
    else:
        if isinstance(lv, str):
            if re.match(level_rgx_str, lv):
                use_level = lv
            else:
                use_level = "MAX"
        else:
            use_level = "MAX"
    return use_level


def get_style_level(
    export_style: StyleType, export_level: Union[str, LvType]
) -> (str, str):
    to_level = get_level(export_level)
    if export_style == StyleType.lipidlynxx:
        export_style = "LipidLynxX"
    elif export_style == StyleType.comp_db:
        export_style = "COMP_DB"
        to_level = "B1"
    else:
        export_style = "LipidLynxX"

    return export_style, to_level


def get_url_safe_str(data: dict) -> str:
    data_json: str = json.dumps(data)
    data_bytes: bytes = base64.urlsafe_b64encode(data_json.encode("utf-8"))
    data_str: str = data_bytes.decode("utf-8")

    return data_str


if __name__ == "__main__":
    x = get_level("B1")
    print(x)
