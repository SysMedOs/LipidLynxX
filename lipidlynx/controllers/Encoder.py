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

from lipidlynx.models.defaults import (
    lipid_class_alias_info,
    cv_order_list,
    cv_alias_info,
)

# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from lipidlynx.models.log import logger
from lipidlynx.controllers.general_functions import seg_to_str
from lipidlynx.controllers.Parser import parse, parse_mod


def lynx_encode(parsed_dct: Dict[str, Union[str, dict, list]]) -> str:
    """
    The general encoder to generate lion code from parsed information
    Args:
        parsed_dct:

    Returns:

    """
    input_abbr = parsed_dct.get("INPUT_ABBR", None)
    logger.debug(f"Try to encode lipid: {input_abbr}")
    lynx_code_candidates_lst = []
    parsed_dct["CLASS_INFO"] = []
    for reg_pattern in parsed_dct["OUTPUT_INFO"]:
        tmp_parsed_dct = parsed_dct["OUTPUT_INFO"][reg_pattern]
        lipid_class = tmp_parsed_dct.get("RULE_CLASS", "")

        if lipid_class is not None:
            if lipid_class in "FA":
                lynx_code_candidates_lst.append(encode_fa(tmp_parsed_dct))
            elif lipid_class == "SPB":
                lynx_code_candidates_lst.append(encode_spb(tmp_parsed_dct))
            elif lipid_class == "PL":
                lynx_code_candidates_lst.append(encode_pl(tmp_parsed_dct))
            elif lipid_class in ["Cer", "SM", "SP"]:
                lynx_code_candidates_lst.append(encode_sp(tmp_parsed_dct))
            elif lipid_class == "GL":
                lynx_code_candidates_lst.append(encode_gl(tmp_parsed_dct))
            elif lipid_class == "CL":
                lynx_code_candidates_lst.append(encode_cl(tmp_parsed_dct))
            elif lipid_class == "BMP":
                lynx_code_candidates_lst.append(encode_bmp(tmp_parsed_dct))

        parsed_dct["CLASS_INFO"].append(tmp_parsed_dct.get("PARSED_CLASS", ""))

    class_info_lst = parsed_dct["CLASS_INFO"]
    try:
        best_class = str(Counter(class_info_lst).most_common(1)[0][0])
        logger.debug(f"Best class is: {best_class}")
    except IndexError:
        best_class = None
        logger.warning(f"Failed to propose best class for this lipid.")
    lynx_code = get_best_abbreviation(
        lynx_code_candidates_lst, target_class=best_class, abbr=input_abbr
    )

    add_mod = parsed_dct.get("ADDITIONAL_MOD", None)
    if add_mod:
        add_mod_lst = decode_mod(add_mod)
        if add_mod_lst:
            add_mod_str = seg_to_str(add_mod_lst).strip("[]")
            add_mod_str = seg_to_str(add_mod_lst).strip("<>")
            lynx_code += f"<{add_mod_str}>"

    return lynx_code


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

    sum_lynx_code_dct = {}

    for p in sum_parsed_info:

        if isinstance(sum_parsed_info, list) or isinstance(sum_parsed_info, tuple):
            parsed_dct = p
        elif isinstance(sum_parsed_info, dict):
            parsed_dct = sum_parsed_info[p]
        else:
            raise TypeError(f"Can NOT process input type {type(sum_parsed_info)}")
        lynx_code = lynx_encode(parsed_dct)
        abbr = parsed_dct.get("INPUT_ABBR", None)
        if not abbr:
            abbr = lynx_code

        sum_lynx_code_dct[abbr] = lynx_code

    return sum_lynx_code_dct


def get_best_abbreviation(
    candidates_lst: List[str], target_class: str = None, abbr: str = None
) -> str:
    lynx_code = ""
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
                    lynx_code = tmp_lion
            else:
                if tmp_len > code_len:
                    code_len = tmp_len
                    lynx_code = tmp_lion
                elif tmp_len == code_len and tmp_lion.startswith(("O-", "P-")):
                    code_len = tmp_len
                    lynx_code = tmp_lion

    elif len(candidates_lst) == 1:
        lynx_code = candidates_lst[0]
    else:
        logger.warning("Failed to generate epiLION abbreviation for this lipid...")

    if lynx_code is None:
        lynx_code = ""

    return lynx_code


def decode_mod(mod_info: str, front: str = "position", end: str = "count") -> List[str]:
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
        cv_replace_lst = []
        cv_count = 0

        for cv_info in cv_info_lst:
            tmp_front = cv_info.get("FRONT", None)
            tmp_end = cv_info.get("END", None)
            tmp_replace = cv_info.get("REPLACE", None)
            if tmp_front and re.match(r"\d\d?[xX]$", tmp_front):
                tmp_front = tmp_front[:-1]
                front = "count"
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

            if tmp_replace:
                rep_pos = re.match(r"@[CHNOP](?P<REPLACE>\d\d?)", tmp_replace)
                if rep_pos:
                    cv_replace_lst.append(rep_pos.groupdict()["REPLACE"])

        if cv_count == 1:
            cv_count_str = ""
        else:
            cv_count_str = str(cv_count)
        if cv_position_lst:
            position_str = "{" + seg_to_str(cv_position_lst) + "}"
            mod_encode_dct[cv] = f"{cv_count_str}{cv}{position_str}"
        elif cv_replace_lst:
            replace_str = "@{" + seg_to_str(cv_replace_lst) + "}"
            mod_encode_dct[cv] = f"{cv_count_str}{cv}{replace_str}"
        else:
            mod_encode_dct[cv] = f"{cv_count_str}{cv}"

    # Sort modifications by order
    encoded_mod_lst = []

    for c in cv_order_list:
        if c in mod_encode_dct:
            encoded_mod_lst.append(mod_encode_dct[c])
    # encoded_mod_str = seg_to_str(encoded_mod_lst)

    return encoded_mod_lst


def get_mod_code(parsed_info: dict, add_mod: str = None, brackets: bool = True) -> str:

    mod_code = ""
    mod_lst = []
    c_count = parsed_info.get("C", None)
    db_info = parsed_info.get("DB_INFO", None)
    mod_info = parsed_info.get("MOD_INFO", None)
    o_info = parsed_info.get("O_INFO", None)
    mod_info_lst = []

    if db_info is not None:
        mod_lst.append("DB{" + db_info.strip("()[]<>") + "}")
    if mod_info is not None:
        mod_info_lst.append(mod_info)
    if add_mod is not None:
        mod_info_lst.append(add_mod)
    if mod_info_lst:
        for mod_str in mod_info_lst:
            if mod_str.strip("()") in ["COOH", "CHO"]:
                rep_str = "@{" + c_count.strip("()[]<>") + "}"
                seg_mod_lst = decode_mod(mod_str)
                seg_code_lst = []
                for seg_mod in seg_mod_lst:
                    if "@{" in seg_mod:
                        seg_code_lst.append(seg_mod)
                    else:
                        seg_code_lst.append(seg_mod + rep_str)
                mod_lst.extend(seg_code_lst)
            else:
                mod_lst.extend(decode_mod(mod_str))

    if o_info is not None:
        o_info_len = len(o_info)
        if o_info_len <= 2 and re.match(r"\d{1,2}", o_info):
            mod_lst.append(f"{o_info}O")
        elif 2 <= o_info_len <= 3 and re.match(r"\d{1,2}O", o_info):
            mod_lst.append(o_info)
        elif o_info_len > 3 and re.match(r"\d{1,2}\(.*\)", o_info):
            mod_lst.extend(decode_mod(o_info[1:].strip("()[]<>")))
        else:
            mod_lst.extend(decode_mod(o_info))
    sorted_mod_lst = []
    mod_lst = [m for m in mod_lst if m is not None]
    mod_lst = [m for m in mod_lst if m != ""]
    if mod_lst:
        for mod in cv_order_list:
            if mod in cv_alias_info:
                mod_alias_lst = cv_alias_info[mod]
                for mod_alia in mod_alias_lst:
                    for obs_mod in mod_lst:
                        if obs_mod and re.match(
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
            mod_code = f"<{mod_code}>"

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

    lynx_code = f'{class_code}{parsed_info["C"]}:{parsed_info["DB"]}'
    lynx_code += get_mod_code(parsed_info, add_mod)
    parsed_info["PARSED_CLASS"] = class_code

    return lynx_code


def encode_spb(parsed_info: dict, add_mod: str = None, is_sub: bool = False) -> str:
    logger.debug(parsed_info)
    link_prefix = parsed_info.get("CLASS", None)
    if link_prefix:
        link_prefix = link_prefix.upper()
        if link_prefix == "SPBP":
            link_prefix = "SPBP<PO4>"
    else:
        link_prefix = "SPB"

    lynx_code = f'{parsed_info["C"]}:{parsed_info["DB"]}'
    lynx_code += get_mod_code(parsed_info, add_mod)

    if is_sub:
        pass
    else:
        lynx_code = f"{link_prefix}({lynx_code})"

    return lynx_code


def encode_sub_fa(fa_abbr: str, add_mod: str = None):
    fa_candidates_lst = []
    class_info_lst = []
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
    try:
        best_class = str(Counter(class_info_lst).most_common(1)[0][0])
    except IndexError:
        best_class = ""
    logger.debug(f"Best class for {fa_abbr} is: {best_class}")

    fa_lynx_code = get_best_abbreviation(
        fa_candidates_lst, target_class=best_class, abbr=fa_abbr
    )
    fa_lynx_code = fa_lynx_code.strip("FA")

    return fa_lynx_code


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
            mod_op_lst = []
            unmod_op_lst = []
            mod_fa_lst = []
            unmod_fa_lst = []

            for fa in fa_info_lst:
                if fa.startswith(("O-", "P-")):
                    if "[" in fa:
                        mod_op_lst.append(fa)
                    elif "<" in fa:
                        mod_op_lst.append(fa)
                    else:
                        unmod_op_lst.append(fa)
                elif "[" in fa:
                    mod_fa_lst.append(fa)
                elif "<" in fa:
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
    lynx_code = f"{lyso_prefix}{class_code}{fa_code}"

    if len(check_empty_fa(lynx_code)) == 1 and not lynx_code.startswith("L"):
        lynx_code = f"L{lynx_code}"

    return lynx_code


def encode_sp(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    spb_abbr = parsed_info["SPB"]
    if not spb_abbr.startswith("SPB"):
        spb_abbr = f"SPB{spb_abbr}"
    spb_code = encode_sub_fa(spb_abbr, add_mod=parsed_info.get("SPB_MOD", None))
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=1, brackets=False)
    class_code = check_lipid_class_alias(parsed_info["CLASS"])
    if fa_code:
        lynx_code = f"{class_code}({spb_code}/{fa_code})"
    else:
        lynx_code = f"{class_code}({spb_code})"
    return lynx_code


def encode_gl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    fa_code = encode_all_sub_fa(parsed_info=parsed_info, fa_count=3)
    class_code = check_lipid_class_alias(parsed_info["CLASS"])
    lynx_code = f"{class_code}{fa_code}"

    if len(check_empty_fa(lynx_code)) == 1 and not lynx_code.startswith("DG"):
        lynx_code = f"DG{lynx_code[3:]}"
    elif len(check_empty_fa(lynx_code)) == 2 and not lynx_code.startswith("MG"):
        lynx_code = f"MG{lynx_code[3:]}"

    return lynx_code


def encode_cl(parsed_info: dict) -> str:
    logger.debug(parsed_info)
    lynx_code = parsed_info["CLASS"] + encode_all_sub_fa(
        parsed_info=parsed_info, fa_count=4
    )

    return lynx_code


def encode_bmp(parsed_info: dict) -> str:
    return encode_cl(parsed_info)


if __name__ == "__main__":

    # x = encode_mod("+2O")
    # print(x)
    # y = encode_mod("+O2")
    # print(y)
    # z = lynx_encode(parse("PE O-18:1p/18:1"))
    z = lynx_encode(parse(r"FA (18:2) + O"))
    print(z)
