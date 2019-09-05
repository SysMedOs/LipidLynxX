# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List, Union

from epilion.controllers.Logger import logger
from epilion.controllers.DefaultParams import class_rgx_dct, rgx_class_dct, cv_rgx_dct


def parse(
    abbr: str,
    class_rules_dct: Dict[str, re.compile] = class_rgx_dct,
    rules_class_dct: Dict[re.compile, str] = rgx_class_dct,
    lipid_class: str = None,
) -> Dict[str, Union[str, dict]]:

    """
    Main parser to read input abbreviations
    Args:
        abbr: input lipid abbreviation to be converted
        class_rules_dct: the predefined dict in the form of lipid_class: re.compile(r"rule")
        rules_class_dct: the predefined dict in the form of re.compile(r"rule"): lipid_class
        lipid_class: name of lipid_class

    Returns:
        parsed_info_dct: parsed information stored as dict

    """

    rgx_lst = []
    abbr = abbr.strip('"')  # remove possible side " from csv
    rgx_white = re.compile(r"\s+")
    abbr = re.sub(rgx_white, "", abbr)  # remove all white space

    if not lipid_class:  # try to get lipid class from abbr
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
            if not rgx_lst:
                rgx_lst = class_rules_dct.get("GL", [])
        elif abbr.upper().startswith(("DG", "DAG")):
            rgx_lst = class_rules_dct.get("DG", [])
            if not rgx_lst:
                rgx_lst = class_rules_dct.get("GL", [])
        elif abbr.upper().startswith(("TG", "TAG")):
            rgx_lst = class_rules_dct.get("TG", [])
            if not rgx_lst:
                rgx_lst = class_rules_dct.get("GL", [])
    else:
        if lipid_class in class_rules_dct:
            rgx_lst = class_rules_dct.get(lipid_class, [])
        else:
            raise ValueError(f"Lipid class {lipid_class} is not supported")

    if not rgx_lst:
        rgx_lst = [r for k in class_rules_dct for r in class_rules_dct[k]]

    parsed_info_dct = get_match_info(abbr, rgx_lst, rules_class_dct=rules_class_dct)

    if not parsed_info_dct:
        logger.warning(
            f'Can not parse abbreviation: "{abbr}", try to ignore case and try again...'
        )
        parsed_info_dct = get_match_info(abbr, rgx_lst, rules_class_dct=rules_class_dct)
        if parsed_info_dct:
            logger.info(f'Successfully parsed abbreviation: "{abbr}"')
        else:
            logger.warning(f'Not able to parse abbreviation: "{abbr}"')

    return parsed_info_dct


def get_match_info(
    abbr: str,
    rgx_lst: List[re.compile],
    rules_class_dct: Dict[re.compile, str] = rgx_class_dct,
    ignore_case: bool = False,
) -> Dict[str, Union[str, dict]]:
    """
    General match function to find pattern by regular expression
    Args:
        abbr: input lipid abbreviation to be converted
        rgx_lst: list of possible regular expression patterns for the input abbr in the form of List[re.compile]
        rules_class_dct: the predefined dict in the form of re.compile(r"rule"): lipid_class
        ignore_case: set to False by default to be strict with cases defined in the patterns

    Returns:
        the matched information as dict

    """
    matched_info_dct = {"INPUT_ABBR": abbr, "OUTPUT_INFO": {}}
    if rgx_lst:
        for rgx in rgx_lst:
            if not ignore_case:
                pass
            else:
                rgx = re.compile(rgx.pattern, flags=re.IGNORECASE)
            if rgx.match(abbr):
                rgx_match = rgx.match(abbr)
                rgx_pattern = str(rgx.pattern)
                matched_info_dct["OUTPUT_INFO"][rgx_pattern] = rgx_match.groupdict()
                if (
                    matched_info_dct["OUTPUT_INFO"][rgx_pattern]
                    and matched_info_dct["OUTPUT_INFO"][rgx_pattern].get("CLASS", None)
                    is None
                ):
                    if rgx in rules_class_dct:
                        matched_info_dct["OUTPUT_INFO"][rgx_pattern][
                            "CLASS"
                        ] = rules_class_dct[rgx]
                    else:
                        raise ValueError(f"Can not get the lipid class of {abbr}")
                else:
                    del matched_info_dct["OUTPUT_INFO"][rgx_pattern]

    return matched_info_dct


def parse_mod(
    abbr: str, cv_patterns_dct: Dict[str, re.compile] = cv_rgx_dct
) -> Dict[str, list]:
    """
    parse the modifications based on predefined list of abbreviations
    Args:
        abbr: input lipid abbreviation to be converted
        cv_patterns_dct: the predefined dict in the form of CV: re.compile(r"rule")

    Returns:
        the matched information as dict

    """
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
