# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List, Tuple, Union

from epilion.controllers.DefaultParams import cv_lst
from epilion.controllers.Logger import logger
from epilion.controllers.GeneralFunctions import seg_to_str
from epilion.controllers.Parser import parse, parse_mod


def lion_encode(parsed_dct: Dict[str, Union[str, dict]]) -> str:
    """
    The general encoder to generate lion code from parsed information
    Args:
        parsed_dct:

    Returns:

    """

    logger.debug(f'Try to encode lipid: {parsed_dct.get("INPUT_ABBR", None)}')
    lion_code_candidates_lst = []
    for reg_pattern in parsed_dct["OUTPUT_INFO"]:
        tmp_parsed_dct = parsed_dct["OUTPUT_INFO"][reg_pattern]
        lipid_class = tmp_parsed_dct.get("CLASS", None)

        if lipid_class is not None:
            if lipid_class == "FA":
                lion_code_candidates_lst.append(encode_fa(tmp_parsed_dct))
            elif lipid_class == "PL":
                lion_code_candidates_lst.append(encode_pl(tmp_parsed_dct))
            elif lipid_class == "GL":
                lion_code_candidates_lst.append(encode_gl(tmp_parsed_dct))
            elif lipid_class == "CL":
                lion_code_candidates_lst.append(encode_cl(tmp_parsed_dct))
            elif lipid_class == "BMP":
                lion_code_candidates_lst.append(encode_bmp(tmp_parsed_dct))

    lion_code = get_best_abbreviation(lion_code_candidates_lst)

    return lion_code


def batch_encode(
    sum_parsed_info: Union[
        Dict[any, Dict[str, Union[str, dict]]],
        List[Dict[str, Union[str, dict]]],
        Tuple[Dict[str, Union[str, dict]]],
    ]
) -> Dict[str, str]:
    """
    Load sum of parsed info to generate a dictionary of lion code output
    Args:
        sum_parsed_info:

    Returns:

    """

    sum_lion_code_dct = {}

    for p in sum_parsed_info:

        if isinstance(sum_parsed_info, list) or isinstance(sum_parsed_info, tuple):
            parsed_dct = p
        elif isinstance(sum_parsed_info, dict):
            parsed_dct = sum_parsed_info[p]
        else:
            parsed_dct = {}
            raise TypeError(f"Can NOT process input type {type(sum_parsed_info)}")
        lion_code = lion_encode(parsed_dct)
        abbr = parsed_dct.get("INPUT_ABBR", None)
        if not abbr:
            abbr = lion_code

        sum_lion_code_dct[abbr] = lion_code

    return sum_lion_code_dct


def get_best_abbreviation(candidates_lst: List[str]) -> str:
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

    if lion_code is None:
        lion_code = ""

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
                if cv_info["END"] and end.lower() == "count":
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


def encode_fa(parsed_info: dict, add_mod: str = None) -> str:
    logger.debug(parsed_info)
    link_prefix = parsed_info.get("LINK", None)
    if link_prefix is not None:
        link_prefix = link_prefix.upper()
    else:
        link_prefix = ""
    lion_code = f'{link_prefix}{parsed_info["C"]}:{parsed_info["DB"]}'
    db_info = parsed_info.get("DB_INFO", None)
    mod_info = parsed_info.get("MOD_INFO", None)
    o_info = parsed_info.get("O_INFO", None)
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
    if fa_abbr:
        fa_dct = parse(fa_abbr)["OUTPUT_INFO"]
        for reg_pattern in fa_dct:
            tmp_parsed_dct = fa_dct[reg_pattern]
            lipid_class = tmp_parsed_dct.get("CLASS", "")
            if lipid_class == "FA":
                fa_candidates_lst.append(encode_fa(tmp_parsed_dct, add_mod=add_mod))
    fa_lion_code = get_best_abbreviation(fa_candidates_lst).strip("FA")

    return fa_lion_code


def encode_all_sub_fa(parsed_info: dict, fa_count: int, brackets: bool = True) -> str:

    fa_count_lst = list(range(1, fa_count + 1))
    position_lst = list(range(1, fa_count))

    fa_info_lst = []

    for f in fa_count_lst:
        fa_abbr = parsed_info.get(f"FA{f}", "")
        if fa_abbr:
            fa_code = encode_sub_fa(fa_abbr, add_mod=parsed_info.get(f"FA{f}_MOD", ""))

        else:
            fa_code = ""
        fa_info_lst.append(fa_code)
        if f in position_lst:
            fa_info_lst.append(parsed_info.get(f"POSITION{f}", ""))
    if brackets:
        fa_code = f'({seg_to_str(fa_info_lst,sep="")})'
    else:
        fa_code = f'{seg_to_str(fa_info_lst, sep="")}'

    return fa_code


def check_fa(abbr: str) -> list:
    no_fa_rgx = re.compile(r"[_/(\-]0:0")
    no_fa_lst = no_fa_rgx.findall(abbr)
    return no_fa_lst


def encode_pl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    lyso_prefix = parsed_info.get("LYSO", None)
    if lyso_prefix is not None:
        lyso_prefix = lyso_prefix.upper()
    else:
        lyso_prefix = ""
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=2)
    lion_code = f'{lyso_prefix}{parsed_info["PL"].upper()}{fa_code}'

    if len(check_fa(lion_code)) == 1 and not lion_code.startswith("L"):
        lion_code = f"L{lion_code}"

    return lion_code


def encode_gl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=3)
    lion_code = f'{parsed_info["GL"].upper()}{fa_code}'

    if len(check_fa(lion_code)) == 1 and not lion_code.startswith("DG"):
        lion_code = f"DG{lion_code[3:]}"
    elif len(check_fa(lion_code)) == 2 and not lion_code.startswith("MG"):
        lion_code = f"MG{lion_code[3:]}"

    return lion_code


def encode_cl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    lion_code = "CL" + encode_all_sub_fa(parsed_info=parsed_info, fa_count=4)

    return lion_code


def encode_bmp(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    lion_code = "BMP" + encode_all_sub_fa(parsed_info=parsed_info, fa_count=4)
    return lion_code
