# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re

from epilion.controllers.DefaultParams import cv_lst
from epilion.controllers.Parser import parse, parse_mod


def lion_encode(parsed_dct: dict) -> str:

    lion_code = ""
    lipid_class = parsed_dct.get("CLASS", None)
    if lipid_class is not None:
        if lipid_class == "FA":
            lion_code = encode_fa(parsed_dct)

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

        if cv_position_lst:
            position_str = "{" + ",".join(cv_position_lst) + "}"
            mod_encode_dct[cv] = f"{cv_count}{cv}{position_str}"
        else:
            mod_encode_dct[cv] = f"{cv_count}{cv}"

    if "O" in mod_encode_dct:
        for o_mod in ["OH", "oxo", "oxo", "ep", "OO", "OOH"]:
            if o_mod in mod_encode_dct:
                mod_encode_dct.pop("O", None)
    encoded_mod_lst = []
    for c in cv_lst:
        if c in mod_encode_dct:
            encoded_mod_lst.append(mod_encode_dct[c])
    encoded_mod_str = ",".join(encoded_mod_lst)

    return encoded_mod_str


def encode_fa(parsed_dct: dict) -> str:
    lion_code = f'{parsed_dct["LINK"].upper()}{parsed_dct["C"]}:{parsed_dct["DB"]}'
    db_info = parsed_dct.get("DB_INFO", None)
    mod_info = parsed_dct.get("MOD_INFO", None)
    o_info = parsed_dct.get("O_INFO", None)
    mod_lst = []
    if db_info is not None:
        mod_lst.append("{" + db_info.strip("()[]") + "}")
    if mod_info is not None:
        mod_lst.append(encode_mod(mod_info))
    if o_info is not None:
        mod_lst.append(encode_mod(o_info))
    if mod_lst:
        lion_code += "[" + ",".join(mod_lst) + "]"

    return lion_code


if __name__ == "__main__":
    usr_rules_file = r"../configurations/rules.csv"

    usr_abbr_lst = [
        r"FA 20:4;O2",
        r"fa 20:4;O2",
        r"Test",
        "FA 20:4_OH2",
        "FA 20:4(6E,8E,10E,14E)(5OH,12OH)",
    ]

    for usr_abbr in usr_abbr_lst:
        usr_parsed_info_dct = parse(usr_abbr)
        print(usr_parsed_info_dct)
        usr_lion_code = lion_encode(usr_parsed_info_dct)
        print(usr_abbr, " -> ", usr_lion_code)
