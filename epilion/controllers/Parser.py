# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re

from epilion.controllers.Logger import logger
from epilion.controllers.DefaultParams import class_rgx_dct, rgx_class_dct


def parse(
    abbr: str,
    class_rules_dct: dict = class_rgx_dct,
    rules_class_dct: dict = rgx_class_dct,
    lipid_class: str = None,
) -> dict:

    parsed_info_dct = {}
    rgx_lst = []
    detected_lipid_class = None

    if not lipid_class:
        if abbr.upper().startswith(("FA", "O-", "P-")):
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

    if rgx_lst:
        for rgx in rgx_lst:
            if rgx.match(abbr):
                rgx_match = rgx.match(abbr)
                parsed_info_dct = rgx_match.groupdict()
                if parsed_info_dct and parsed_info_dct.get("CLASS", None) is None:
                    if rgx in rules_class_dct:
                        parsed_info_dct["CLASS"] = rules_class_dct[rgx]

    if not parsed_info_dct:
        logger.warning(
            f'Can not parse abbreviation: "{abbr}", try to ignore case and try again...'
        )
        for rgx in rgx_lst:
            rgx_no_case = re.compile(rgx.pattern, flags=re.IGNORECASE)
            if rgx_no_case.match(abbr):
                rgx_match = rgx_no_case.match(abbr)
                parsed_info_dct = rgx_match.groupdict()
                if parsed_info_dct and parsed_info_dct.get("CLASS", None) is None:
                    if rgx in rules_class_dct:
                        parsed_info_dct["CLASS"] = rules_class_dct[rgx]
        if parsed_info_dct:
            logger.info(f'Successfully parsed abbreviation: "{abbr}"')
        else:
            logger.warning(f'Notable to parse abbreviation: "{abbr}"')

    return parsed_info_dct


if __name__ == "__main__":

    usr_abbr_lst = [r"FA 20:4;O2", r"fa 20:4;O2", r"Test"]

    for usr_abbr in usr_abbr_lst:
        usr_parsed_info_dct = parse(usr_abbr)
        print(usr_parsed_info_dct)
