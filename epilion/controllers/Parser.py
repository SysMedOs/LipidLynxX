# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import List

from epilion.controllers.Logger import logger
from epilion.controllers.DefaultParams import class_rgx_dct, rgx_class_dct, cv_rgx_dct


def parse(
    abbr: str,
    class_rules_dct: dict = class_rgx_dct,
    rules_class_dct: dict = rgx_class_dct,
    lipid_class: str = None,
) -> dict:

    parsed_info_dct = {}
    rgx_lst = []
    detected_lipid_class = None
    abbr = abbr.strip('"')  # remove possible side " from csv
    rgx_white = re.compile(r"\s+")
    abbr = re.sub(rgx_white, "", abbr)  # remove all white space

    if not lipid_class:
        if abbr.upper().startswith(("FA", "O-", "P-")):
            rgx_lst = class_rules_dct.get("FA", [])
        elif abbr[1:].upper().startswith(("O-", "P-")):  # for dO-, dP-
            rgx_lst = class_rules_dct.get("FA", [])
        elif abbr.upper().startswith(("PL", "PA", "PC", "PE", "PG", "PS", "PI")):
            rgx_lst = class_rules_dct.get("PL", [])
        elif abbr.upper().startswith(("SM", "SP")):
            rgx_lst = class_rules_dct.get("SM", [])
        elif abbr.upper().startswith(("CER")):
            rgx_lst = class_rules_dct.get("Cer", [])
        elif abbr.upper().startswith(("MG", "MAG")):
            rgx_lst = class_rules_dct.get("MG", [])
        elif abbr.upper().startswith(("DG", "DAG")):
            rgx_lst = class_rules_dct.get("DG", [])
        elif abbr.upper().startswith(("TG", "TAG")):
            rgx_lst = class_rules_dct.get("TG", [])
    else:
        if lipid_class in class_rules_dct:
            rgx_lst = class_rules_dct.get(lipid_class, [])

    if not rgx_lst:
        rgx_lst = [r for k in class_rules_dct for r in class_rules_dct[k]]

    parsed_info_dct = get_match_info(
        abbr, rgx_lst, parsed_info_dct, rules_class_dct=rules_class_dct
    )

    if not parsed_info_dct:
        logger.warning(
            f'Can not parse abbreviation: "{abbr}", try to ignore case and try again...'
        )
        parsed_info_dct = get_match_info(
            abbr, rgx_lst, parsed_info_dct, rules_class_dct=rules_class_dct
        )
        if parsed_info_dct:
            logger.info(f'Successfully parsed abbreviation: "{abbr}"')
        else:
            logger.warning(f'Not able to parse abbreviation: "{abbr}"')

    return parsed_info_dct


def get_match_info(
    abbr: str,
    rgx_lst: List[re.compile],
    parsed_info_dct: dict,
    rules_class_dct: dict = rgx_class_dct,
    ignore_case: bool = False,
) -> dict:

    if rgx_lst:

        for rgx in rgx_lst:
            if not ignore_case:
                pass
            else:
                rgx = re.compile(rgx.pattern, flags=re.IGNORECASE)
            if rgx.match(abbr):
                rgx_match = rgx.match(abbr)
                parsed_info_dct[rgx.pattern] = rgx_match.groupdict()
                if (
                    parsed_info_dct[rgx.pattern]
                    and parsed_info_dct[rgx.pattern].get("CLASS", None) is None
                ):
                    if rgx in rules_class_dct:
                        parsed_info_dct[rgx.pattern]["CLASS"] = rules_class_dct[rgx]

    return parsed_info_dct


def parse_mod(abbr: str, cv_patterns_dct: dict = cv_rgx_dct) -> dict:
    mod_dct = {}
    for cv in cv_patterns_dct:
        rgx = cv_patterns_dct[cv]
        found_segments_lst = rgx.findall(abbr)
        if found_segments_lst:
            for segment_tpl in found_segments_lst:
                segment_str = "".join(segment_tpl)
                m = rgx.search(segment_str)
                if m:
                    g = m.groupdict()
                    if cv not in mod_dct:
                        mod_dct[g["MOD"]] = [g]
                    else:
                        mod_dct[g["MOD"]].append(g)

    return mod_dct


if __name__ == "__main__":

    usr_abbr_lst = [
        r"FA 20:4;O2",
        r"fa 20:4;O2",
        r"Test",
        "FA 20:4(6E,8E,10E,14E)(5OH,12OH)",
    ]

    for usr_abbr in usr_abbr_lst:
        usr_parsed_info_dct = parse(usr_abbr)
        print(usr_parsed_info_dct)
