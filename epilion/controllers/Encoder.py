# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List, Tuple, Union

from natsort import natsorted

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
        lipid_class = tmp_parsed_dct.get("RULE_CLASS", "")

        if lipid_class is not None:
            if lipid_class in "FA":
                lion_code_candidates_lst.append(encode_fa(tmp_parsed_dct))
            elif lipid_class == "SPB":
                lion_code_candidates_lst.append(encode_spb(tmp_parsed_dct))
            elif lipid_class == "PL":
                lion_code_candidates_lst.append(encode_pl(tmp_parsed_dct))
            elif lipid_class in ["Cer", "SM", "SP"]:
                lion_code_candidates_lst.append(encode_sp(tmp_parsed_dct))
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
    link_prefix = parsed_info.get("CLASS", None)
    if link_prefix:
        link_prefix = link_prefix.upper()
    else:
        link_prefix = "FA"

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
        if re.match(r"\d{1,2}", o_info):
            mod_lst.append(f"{o_info}O")
        else:
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


def encode_spb(parsed_info: dict, add_mod: str = None, is_sub: bool = False) -> str:
    logger.debug(parsed_info)
    link_prefix = parsed_info.get("CLASS", None)
    if link_prefix:
        link_prefix = link_prefix.upper()
        if link_prefix == "SPBP":
            link_prefix = "SPBP[PO4]"
    else:
        link_prefix = "SPB"

    lion_code = f'{parsed_info["C"]}:{parsed_info["DB"]}'
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
        o_info_len = len(o_info)
        if o_info_len <= 2 and re.match(r"\d{1,2}", o_info):
            mod_lst.append(f"{o_info}O")
        elif 2 <= o_info_len <= 3 and re.match(r"\d{1,2}O", o_info):
            mod_lst.append(o_info)
        elif o_info_len > 3 and re.match(r"\d{1,2}\(.*\)", o_info):
            mod_lst.append(encode_mod(o_info[1:].strip("()[]")))
        else:
            mod_lst.append(encode_mod(o_info))
    if mod_lst:
        mod_str = seg_to_str(mod_lst)
        if mod_str:
            lion_code += f"[{mod_str}]"

    if is_sub:
        pass
    else:
        lion_code = f"{link_prefix}({lion_code})"

    return lion_code


def encode_sub_fa(fa_abbr: str, add_mod: str = None):
    fa_candidates_lst = []
    if fa_abbr:
        fa_dct = parse(fa_abbr)["OUTPUT_INFO"]
        for reg_pattern in fa_dct:
            tmp_parsed_dct = fa_dct[reg_pattern]
            lipid_class = tmp_parsed_dct.get("RULE_CLASS", "")
            if lipid_class == "FA":
                fa_candidates_lst.append(encode_fa(tmp_parsed_dct, add_mod=add_mod))
            if lipid_class == "SPB":
                fa_candidates_lst.append(
                    encode_spb(tmp_parsed_dct, add_mod=add_mod, is_sub=True)
                )
    fa_lion_code = get_best_abbreviation(fa_candidates_lst).strip("FA")

    return fa_lion_code


def encode_all_sub_fa(
    parsed_info: dict, fa_count: int, brackets: bool = True, sort_discrete: bool = True
) -> str:

    fa_count_lst = list(range(1, fa_count + 1))
    position_lst = list(range(1, fa_count))

    fa_info_lst = []
    position_info_lst = []
    encode_info_lst = []

    for f in fa_count_lst:
        fa_abbr = parsed_info.get(f"FA{f}", "")
        if fa_abbr:
            fa_code = encode_sub_fa(fa_abbr, add_mod=parsed_info.get(f"FA{f}_MOD", ""))
        else:
            fa_code = ""
        fa_info_lst.append(fa_code)
        encode_info_lst.append(fa_code)
        if f in position_lst:
            position_info_lst.append(parsed_info.get(f"POSITION{f}", ""))
            encode_info_lst.append(parsed_info.get(f"POSITION{f}", ""))

    if position_info_lst == ["_"] * len(position_info_lst):
        if sort_discrete is True:
            sorted_fa_info_lst = natsorted(fa_info_lst)

            for i in fa_count_lst:
                if sorted_fa_info_lst[-1].startswith(("O-", "P-")):
                    sorted_fa_info_lst = [sorted_fa_info_lst[-1]] + sorted_fa_info_lst[
                        :-1
                    ]
            fa_code = f'{seg_to_str(sorted_fa_info_lst, sep="_")}'
        else:
            fa_code = f'{seg_to_str(fa_info_lst, sep="_")}'
    else:
        fa_code = f'{seg_to_str(encode_info_lst, sep="")}'

    if brackets:
        fa_code = f"({fa_code})"

    return fa_code


def check_fa(abbr: str) -> list:
    no_fa_rgx = re.compile(r"[_/(\-]0:0")
    no_fa_lst = no_fa_rgx.findall(abbr)
    return no_fa_lst


def encode_pl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    lyso_prefix = parsed_info.get("CLASS_PREFIX", None)
    if lyso_prefix:
        lyso_prefix = lyso_prefix.upper()
    else:
        lyso_prefix = ""
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=2)
    lion_code = f'{lyso_prefix}{parsed_info["CLASS"].upper()}{fa_code}'

    if len(check_fa(lion_code)) == 1 and not lion_code.startswith("L"):
        lion_code = f"L{lion_code}"

    return lion_code


def encode_sp(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    spb_abbr = parsed_info["SPB"]
    if not spb_abbr.startswith("SPB"):
        spb_abbr = f"SPB{spb_abbr}"
    spb_code = encode_sub_fa(spb_abbr, add_mod=parsed_info.get("SPB_MOD", None))
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=1, brackets=False)
    class_code = parsed_info["CLASS"].upper()
    if class_code == "CER":
        class_code = "Cer"
    if fa_code:
        lion_code = f"{class_code}({spb_code}/{fa_code})"
    else:
        lion_code = f"{class_code}({spb_code})"
    return lion_code


def encode_gl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=3)
    class_code = parsed_info["CLASS"].upper()
    if len(class_code) == 3 and class_code[1] == "A":
        class_code = class_code[0] + class_code[2]
    lion_code = f"{class_code}{fa_code}"

    if len(check_fa(lion_code)) == 1 and not lion_code.startswith("DG"):
        lion_code = f"DG{lion_code[3:]}"
    elif len(check_fa(lion_code)) == 2 and not lion_code.startswith("MG"):
        lion_code = f"MG{lion_code[3:]}"

    return lion_code


def encode_cl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    lion_code = parsed_info["CLASS"] + encode_all_sub_fa(
        parsed_info=parsed_info, fa_count=4
    )

    return lion_code


def encode_bmp(parsed_info: dict) -> str:
    return encode_cl(parsed_info)
