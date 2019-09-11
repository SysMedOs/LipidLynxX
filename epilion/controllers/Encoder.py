# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from collections import Counter
import re
from typing import Dict, List, Tuple, Union

from natsort import natsorted

from epilion.controllers.DefaultParams import (
    lipid_class_alias_info,
    cv_order_list,
    cv_alias_info,
)
from epilion.controllers.Logger import logger
from epilion.controllers.GeneralFunctions import seg_to_str
from epilion.controllers.Parser import parse, parse_mod


def lion_encode(parsed_dct: Dict[str, Union[str, dict, list]]) -> str:
    """
    The general encoder to generate lion code from parsed information
    Args:
        parsed_dct:

    Returns:

    """
    input_abbr = parsed_dct.get("INPUT_ABBR", None)
    logger.debug(f"Try to encode lipid: {input_abbr}")
    lion_code_candidates_lst = []
    parsed_dct["CLASS_INFO"] = []
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

        parsed_dct["CLASS_INFO"].append(tmp_parsed_dct.get("PARSED_CLASS", ""))

    class_info_lst = parsed_dct["CLASS_INFO"]
    try:
        best_class = str(Counter(class_info_lst).most_common(1)[0][0])
        logger.debug(f"Best class is: {best_class}")
    except IndexError:
        best_class = None
        logger.warning(f"Failed to propose best class for this lipid.")
    lion_code = get_best_abbreviation(
        lion_code_candidates_lst, target_class=best_class, abbr=input_abbr
    )

    add_mod = parsed_dct.get("ADDITIONAL_MOD", None)
    if add_mod:
        add_mod_str = encode_mod(add_mod)
        if add_mod_str:
            lion_code += f"[{add_mod_str}]"

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
            raise TypeError(f"Can NOT process input type {type(sum_parsed_info)}")
        lion_code = lion_encode(parsed_dct)
        abbr = parsed_dct.get("INPUT_ABBR", None)
        if not abbr:
            abbr = lion_code

        sum_lion_code_dct[abbr] = lion_code

    return sum_lion_code_dct


def get_best_abbreviation(
    candidates_lst: List[str], target_class: str = None, abbr: str = None
) -> str:
    lion_code = ""
    candidates_lst = list(filter(None, list(set(candidates_lst))))
    if abbr:
        if re.search(r":\d[oep]", abbr) or re.search(r"[OP]-\d{1,2}:\da?", abbr):
            for c in candidates_lst:
                if c.startswith("FA"):
                    candidates_lst.remove(c)
    if len(candidates_lst) > 1:
        code_len = 1
        for tmp_lion in candidates_lst:
            tmp_len = len(tmp_lion)
            if target_class:
                if tmp_lion.startswith(target_class) and tmp_len > code_len:
                    code_len = tmp_len
                    lion_code = tmp_lion
            else:
                if tmp_len > code_len:
                    code_len = tmp_len
                    lion_code = tmp_lion
                elif tmp_len == code_len and tmp_lion.startswith(("O-", "P-")):
                    code_len = tmp_len
                    lion_code = tmp_lion

    elif len(candidates_lst) == 1:
        lion_code = candidates_lst[0]
    else:
        logger.warning("Failed to generate epiLION abbreviation for this lipid...")

    if lion_code is None:
        lion_code = ""

    return lion_code


def encode_mod(mod_info: str, front: str = "position", end: str = "count") -> str:
    if mod_info:
        mod_dct = parse_mod(mod_info)
    else:
        return ""

    if front:
        if front.lower().startswith("count") or "[" in mod_info:
            front = "count"
        elif mod_info.startswith("+"):
            front = "count"
        else:
            front = "position"
    else:
        front = "position"
    if not end:
        end = "count"

    mod_encode_dct = {}

    for cv in mod_dct:
        cv_info_lst = mod_dct[cv]
        cv_position_lst = []
        cv_count = 0

        for cv_info in cv_info_lst:
            tmp_front = cv_info.get("FRONT", None)
            tmp_end = cv_info.get("END", None)
            if tmp_front and tmp_end:
                cv_position_lst.append(str(tmp_front))
                cv_count += 1
            elif tmp_front and not tmp_end:
                if front == "position":
                    cv_position_lst.append(str(tmp_front))
                    cv_count += 1
                else:
                    try:
                        _count = int(tmp_front)
                    except ValueError:
                        _count = 1
                    cv_count += _count
            elif not tmp_front and tmp_end and end.lower() == "count":
                try:
                    _count = int(tmp_end)
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

    # Sort modifications by order
    encoded_mod_lst = []

    for c in cv_order_list:
        if c in mod_encode_dct:
            encoded_mod_lst.append(mod_encode_dct[c])
    encoded_mod_str = seg_to_str(encoded_mod_lst)

    return encoded_mod_str


def get_mod_code(parsed_info: dict, add_mod: str = None, brackets: bool = True) -> str:

    mod_code = ""
    mod_lst = []
    c_count = parsed_info.get("C", None)
    db_info = parsed_info.get("DB_INFO", None)
    mod_info = parsed_info.get("MOD_INFO", None)
    o_info = parsed_info.get("O_INFO", None)

    if db_info is not None:
        mod_lst.append("DB{" + db_info.strip("()[]") + "}")
    if mod_info is not None:
        if mod_info.strip("()") in ["COOH", "CHO"]:
            mod_lst.append(encode_mod(mod_info) + "@{" + c_count.strip("()[]") + "}")
        else:
            mod_lst.append(encode_mod(mod_info))
    if add_mod is not None:
        if add_mod.strip("()") in ["COOH", "CHO"]:
            mod_lst.append(encode_mod(add_mod) + "@{" + c_count.strip("()[]") + "}")
        else:
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
    sorted_mod_lst = []
    mod_lst = [m for m in mod_lst if m is not None]
    mod_lst = [m for m in mod_lst if m != ""]
    if mod_lst:
        for _obs_mod in mod_lst:
            # check if still any not split modifications
            if "{" not in _obs_mod:
                if "," in _obs_mod or "_" in _obs_mod:
                    if len(_obs_mod.split(",")) > 1:
                        mod_lst.extend(_obs_mod.split(","))
                        mod_lst.remove(_obs_mod)
                    elif len(_obs_mod.split("_")) > 1:
                        mod_lst.extend(_obs_mod.split("_"))
                        mod_lst.remove(_obs_mod)
            elif "}," in _obs_mod:
                try:
                    rep_obs_mod = re.sub(r'},', '}|', _obs_mod)
                    if len(rep_obs_mod.split("|")) > 1:
                        mod_lst.extend(rep_obs_mod.split("|"))
                        mod_lst.remove(_obs_mod)
                except TypeError:
                    pass
        for mod in cv_order_list:
            if mod in cv_alias_info:
                mod_alias_lst = cv_alias_info[mod]
                for mod_alia in mod_alias_lst:
                    for obs_mod in mod_lst:
                        if re.match(
                            r"\d{0,2}%s[@]?({([,]?\d{1,2}[EZez]?)+})?$" % mod_alia,
                            obs_mod,
                        ):
                            if obs_mod.startswith("DB"):
                                obs_mod = obs_mod.strip("DB")
                            sorted_mod_lst.append(obs_mod)
            else:
                logger.warning(f"Modification {mod} do not have defined alias.")

    if sorted_mod_lst:
        mod_code = seg_to_str(sorted_mod_lst)
        if mod_code and brackets:
            mod_code = f"[{mod_code}]"

    return mod_code


def encode_fa(parsed_info: dict, add_mod: str = None) -> str:
    logger.debug(parsed_info)

    # carefully detect if FA is FA , O- or P-
    class_code = parsed_info.get("CLASS", None)
    if class_code:
        class_code = check_fa_link(class_code)
    link = parsed_info.get("LINK", None)
    if link:
        link_code = check_fa_link(link)
    else:
        link_code = None
    if not class_code:
        class_code = "FA"
    if link_code and class_code != link_code:
        class_code = link_code

    lion_code = f'{class_code}{parsed_info["C"]}:{parsed_info["DB"]}'
    lion_code += get_mod_code(parsed_info, add_mod)
    parsed_info["PARSED_CLASS"] = class_code

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
    lion_code += get_mod_code(parsed_info, add_mod)

    if is_sub:
        pass
    else:
        lion_code = f"{link_prefix}({lion_code})"

    return lion_code


def encode_sub_fa(fa_abbr: str, add_mod: str = None):
    fa_candidates_lst = []
    class_info_lst = []
    best_class = None
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
            class_info_lst.append(tmp_parsed_dct.get("PARSED_CLASS", ""))
    best_class = str(Counter(class_info_lst).most_common(1)[0][0])
    logger.debug(f"Best class for {fa_abbr} is: {best_class}")
    fa_lion_code = get_best_abbreviation(
        fa_candidates_lst, target_class=best_class, abbr=fa_abbr
    )
    fa_lion_code = fa_lion_code.strip("FA")

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
            sorted_fa_info_lst = []
            mod_op_lst = []
            unmod_op_lst = []
            mod_fa_lst = []
            unmod_fa_lst = []

            for fa in fa_info_lst:
                if fa.startswith(("O-", "P-")):
                    if "[" in fa:
                        mod_op_lst.append(fa)
                    else:
                        unmod_op_lst.append(fa)
                elif "[" in fa:
                    mod_fa_lst.append(fa)
                else:
                    unmod_fa_lst.append(fa)

            if natsorted(fa_info_lst) == natsorted(
                mod_op_lst + unmod_op_lst + mod_fa_lst + unmod_fa_lst
            ):
                sorted_fa_info_lst = (
                    natsorted(unmod_op_lst)
                    + natsorted(mod_op_lst)
                    + natsorted(unmod_fa_lst)
                    + natsorted(mod_fa_lst)
                )
            else:
                logger.warning("Failed to sort lipid.")
                sorted_fa_info_lst = natsorted(fa_info_lst)

            fa_code = f'{seg_to_str(sorted_fa_info_lst, sep="_")}'
        else:
            fa_code = f'{seg_to_str(fa_info_lst, sep="_")}'
    else:
        fa_code = f'{seg_to_str(encode_info_lst, sep="")}'

    if brackets:
        fa_code = f"({fa_code})"

    return fa_code


def check_empty_fa(abbr: str) -> list:
    no_fa_rgx = re.compile(r"[_/(\-]0:0")
    no_fa_lst = no_fa_rgx.findall(abbr)
    return no_fa_lst


def check_fa_link(link: str) -> str:
    prefix = ""
    if link in ["o", "e", "O-a", "O-e", "O-o", "O-"]:
        prefix = "O-"
    elif link in ["p", "O-p", "P-"]:
        prefix = "P-"
    elif link in ["a"]:
        prefix = "FA"
    elif link.upper().startswith(("MO-", "MP-")):
        prefix = f"{link[1:]}"
    elif link.upper().startswith(("DO-", "DP-")):
        prefix = f"2{link[1:]}"
    elif link.upper().startswith(("TO-", "TP-")):
        prefix = f"3{link[1:]}"
    return prefix


def check_lipid_class_alias(lipid_class: str) -> str:
    if lipid_class in lipid_class_alias_info:
        lipid_class = lipid_class_alias_info[lipid_class]["CLASS"]
    return lipid_class


def encode_pl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    lyso_prefix = parsed_info.get("CLASS_PREFIX", None)
    if lyso_prefix:
        lyso_prefix = lyso_prefix.upper()
    else:
        lyso_prefix = ""
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=2)
    # replace legacy pl class e.g. GPCho to PC
    class_code = check_lipid_class_alias(parsed_info["CLASS"])
    lion_code = f"{lyso_prefix}{class_code}{fa_code}"

    if len(check_empty_fa(lion_code)) == 1 and not lion_code.startswith("L"):
        lion_code = f"L{lion_code}"

    return lion_code


def encode_sp(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    spb_abbr = parsed_info["SPB"]
    if not spb_abbr.startswith("SPB"):
        spb_abbr = f"SPB{spb_abbr}"
    spb_code = encode_sub_fa(spb_abbr, add_mod=parsed_info.get("SPB_MOD", None))
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=1, brackets=False)
    class_code = check_lipid_class_alias(parsed_info["CLASS"])
    if fa_code:
        lion_code = f"{class_code}({spb_code}/{fa_code})"
    else:
        lion_code = f"{class_code}({spb_code})"
    return lion_code


def encode_gl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=3)
    class_code = check_lipid_class_alias(parsed_info["CLASS"])
    lion_code = f"{class_code}{fa_code}"

    if len(check_empty_fa(lion_code)) == 1 and not lion_code.startswith("DG"):
        lion_code = f"DG{lion_code[3:]}"
    elif len(check_empty_fa(lion_code)) == 2 and not lion_code.startswith("MG"):
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


if __name__ == "__main__":

    # x = encode_mod("+2O")
    # print(x)
    # y = encode_mod("+O2")
    # print(y)
    # z = lion_encode(parse("PE O-18:1p/18:1"))
    z = lion_encode(parse("18:2(9Z,12Z)(8OH,11oxo)"))
    print(z)
