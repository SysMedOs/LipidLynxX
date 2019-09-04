# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re

from epilion.controllers.DefaultParams import cv_lst
from epilion.controllers.Logger import logger
from epilion.controllers.GeneralFunctions import seg_to_str
from epilion.controllers.Parser import parse, parse_mod


def lion_encode(parsed_dct: dict) -> str:

    lion_code = ""
    lion_code_candidates_lst = []
    for reg_pattern in parsed_dct:
        tmp_parsed_dct = parsed_dct[reg_pattern]
        lipid_class = tmp_parsed_dct.get("CLASS", None)
        if lipid_class is not None:
            if lipid_class == "FA":
                lion_code_candidates_lst.append(encode_fa(tmp_parsed_dct))

    lion_code_candidates_lst = list(filter(None, list(set(lion_code_candidates_lst))))
    if len(lion_code_candidates_lst) > 1:
        code_len = 1
        for tmp_lion in lion_code_candidates_lst:
            tmp_len = len(tmp_lion)
            if tmp_len > code_len:
                code_len = tmp_len
                lion_code = tmp_lion
    elif len(lion_code_candidates_lst) == 1:
        lion_code = lion_code_candidates_lst[0]
    else:
        logger.warning("Failed to generate epiLION abbreviation for this lipid...")

    return lion_code


def encode_mod(mod_info, front: str = "position", end: str = "count"):
    mod_dct = parse_mod(mod_info)
    if front.lower().startswith("count") or "[" in mod_info:
        front = "count"
    else:
        front = "position"

    mod_encode_dct = {}

    for cv in mod_dct:
        cv_info_lst = mod_dct[cv]
        cv_position_lst = []
        cv_count = 0
        if len(cv_info_lst) == 1 and front == "count":
            cv_info = cv_info_lst[0]
            if cv_info["FRONT"]:
                cv_count += int(cv_info["FRONT"])
            else:
                cv_count += 1
        else:
            for cv_info in cv_info_lst:
                if cv_info["FRONT"]:
                    cv_position_lst.append(str(cv_info["FRONT"]))
                if cv_info["END"]:
                    try:
                        _count = int(cv_info["END"])
                    except ValueError:
                        _count = 1
                    cv_count += _count
                else:
                    cv_count += 1
        if cv_count == 1:
            cv_count_str = ""
        else:
            cv_count_str = str(cv_count)
        if cv_position_lst:
            position_str = "{" + ",".join(cv_position_lst) + "}"
            mod_encode_dct[cv] = f"{cv_count_str}{cv}{position_str}"
        else:
            mod_encode_dct[cv] = f"{cv_count_str}{cv}"

    if "O" in mod_encode_dct:
        for o_mod in ["OH", "oxo", "oxo", "ep", "OO", "OOH"]:
            if o_mod in mod_encode_dct:
                mod_encode_dct.pop("O", None)
    encoded_mod_lst = []
    for c in cv_lst:
        if c in mod_encode_dct:
            encoded_mod_lst.append(mod_encode_dct[c])
    encoded_mod_str = seg_to_str(encoded_mod_lst)

    return encoded_mod_str


def encode_fa(parsed_info_dct: dict) -> str:
    lion_code = f'{parsed_info_dct["LINK"].upper()}{parsed_info_dct["C"]}:{parsed_info_dct["DB"]}'
    db_info = parsed_info_dct.get("DB_INFO", None)
    mod_info = parsed_info_dct.get("MOD_INFO", None)
    o_info = parsed_info_dct.get("O_INFO", None)
    mod_lst = []
    if db_info is not None:
        mod_lst.append("{" + db_info.strip("()[]") + "}")
    if mod_info is not None:
        mod_lst.append(encode_mod(mod_info))
    if o_info is not None:
        mod_lst.append(encode_mod(o_info))
    if mod_lst:
        mod_str = seg_to_str(mod_lst)
        if mod_str:
            lion_code += f"[{mod_str}]"

    return lion_code
