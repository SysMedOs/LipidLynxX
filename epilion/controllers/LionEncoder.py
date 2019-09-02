# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re

from epilion.controllers.Parser import parse


def lion_encode(parsed_dct: dict):

    lion_code = ""
    lipid_class = parsed_dct.get("CLASS", None)
    if lipid_class is not None:
        if lipid_class == "FA":
            lion_code = (
                f'{parsed_dct["LINK"].upper()}{parsed_dct["C"]}:{parsed_dct["DB"]}'
            )
            db_info = parsed_dct.get("DB_INFO", None)
            mod_info = parsed_dct.get("MOD_INFO", None)
            o_info = parsed_dct.get("O_INFO", None)
            mod_lst = []
            if db_info is not None:
                mod_lst.append("{" + db_info.strip("()[]") + "}")
            if mod_info is not None:
                mod_lst.append("{" + mod_info.strip("()[]") + "}")
            if o_info is not None:
                o_rgx = re.compile(r"(?P<O_INFO>[CHON]{0,4})?(?P<O_COUNT>\d)")
                if o_rgx.match(o_info):
                    o_rgx_dict = o_rgx.match(o_info).groupdict()
                    mod_lst.append(f'{o_rgx_dict["O_COUNT"]}{o_rgx_dict["O_INFO"]}')
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
