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
            if lipid_class == "PL":
                lion_code_candidates_lst.append(encode_pl(tmp_parsed_dct))

    lion_code = get_best_abbreviation(lion_code_candidates_lst)

    return lion_code


def get_best_abbreviation(candidates_lst):
    lion_code = ""
    candidates_lst = list(filter(None, list(set(candidates_lst))))
    if len(candidates_lst) > 1:
        code_len = 1
        for tmp_lion in candidates_lst:
            tmp_len = len(tmp_lion)
            if tmp_len > code_len:
                code_len = tmp_len
                lion_code = tmp_lion
    elif len(candidates_lst) == 1:
        lion_code = candidates_lst[0]
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


def encode_fa(parsed_info_dct: dict, add_mod: str = None) -> str:
    link_prefix = parsed_info_dct.get("LINK", None)
    if link_prefix is not None:
        link_prefix = link_prefix.upper()
    else:
        link_prefix = ""
    lion_code = f'{link_prefix}{parsed_info_dct["C"]}:{parsed_info_dct["DB"]}'
    db_info = parsed_info_dct.get("DB_INFO", None)
    mod_info = parsed_info_dct.get("MOD_INFO", None)
    o_info = parsed_info_dct.get("O_INFO", None)
    mod_lst = []
    if db_info is not None:
        mod_lst.append("{" + db_info.strip("()[]") + "}")
    if mod_info is not None:
        mod_lst.append(encode_mod(mod_info))
    if add_mod is not None:
        mod_lst.append(encode_mod(add_mod))
    if o_info is not None:
        mod_lst.append(encode_mod(o_info))
    if mod_lst:
        mod_str = seg_to_str(mod_lst)
        if mod_str:
            lion_code += f"[{mod_str}]"

    if lion_code.upper().startswith(("MO-", "MP-")):
        lion_code = f"{lion_code[1:]}"
    elif lion_code.upper().startswith(("DO-", "DP-")):
        lion_code = f"2{lion_code[1:]}"
    elif lion_code.upper().startswith(("TO-", "TP-")):
        lion_code = f"3{lion_code[1:]}"

    return lion_code


def encode_sub_fa(fa_abbr: str, add_mod: str = None):
    fa_candidates_lst = []
    if fa_abbr is not None:
        fa_dct = parse(fa_abbr)
        for reg_pattern in fa_dct:
            tmp_parsed_dct = fa_dct[reg_pattern]
            lipid_class = tmp_parsed_dct.get("CLASS", None)
            if lipid_class is not None and lipid_class == "FA":
                fa_candidates_lst.append(encode_fa(tmp_parsed_dct, add_mod=add_mod))
    fa_lion_code = get_best_abbreviation(fa_candidates_lst).strip("FA")

    return fa_lion_code


def encode_pl(parsed_info_dct: dict) -> str:
    lyso_prefix = parsed_info_dct.get("LYSO", None)
    if lyso_prefix is not None:
        lyso_prefix = lyso_prefix.upper()
    else:
        lyso_prefix = ""
    lion_code = f'{lyso_prefix}{parsed_info_dct["PL"].upper()}'
    fa1_abbr = parsed_info_dct.get("FA1", None)
    fa2_abbr = parsed_info_dct.get("FA2", None)
    fa1_mod_abbr = parsed_info_dct.get("FA1_MOD", None)
    fa2_mod_abbr = parsed_info_dct.get("FA2_MOD", None)
    position1_info = parsed_info_dct.get("POSITION1", None)
    print(parsed_info_dct)
    fa1_lion_code, fa2_lion_code = "", ""
    if fa1_abbr is not None:
        fa1_lion_code = encode_sub_fa(fa1_abbr, add_mod=fa1_mod_abbr)
    if fa2_abbr is not None:
        fa2_lion_code = encode_sub_fa(fa2_abbr, add_mod=fa2_mod_abbr)
    if position1_info:
        lion_code += f"({fa1_lion_code}{position1_info}{fa2_lion_code})"
    else:
        lion_code += f"({fa1_lion_code})"

    for no_fa in ["/0:0)", "_0:0)", "(0:0_", "(0:0/"]:
        if no_fa in lion_code and not lion_code.startswith("L"):
            lion_code = f"L{lion_code}"

    return lion_code
