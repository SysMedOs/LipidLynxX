# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List, Union

from lipidlynx.models.log import logger
from lipidlynx.models.defaults import (
    class_rgx_dct,
    rgx_class_dct,
    cv_rgx_dct,
    lipid_class_alias_info,
    cv_order_list,
    cv_alias_info,
)
from lipidlynx.controllers.general_functions import seg_to_str


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
    input_abbr = abbr
    abbr = remove_prefix(abbr)

    if "+" in abbr:
        abbr_lst = abbr.split("+")
        abbr = abbr_lst[0]
        abbr_mod = "+" + "+".join(abbr_lst[1:])
        logger.info(f"{input_abbr} contains additional modification {abbr_mod}.")
    else:
        abbr_mod = ""

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
        elif re.match(r"\d{1,2}:.*", abbr):
            rgx_lst = class_rules_dct.get("FA", [])
        elif abbr.startswith(tuple(lipid_class_alias_info.keys())):
            for tmp_class in lipid_class_alias_info:
                if abbr.startswith(tmp_class):
                    rule_class = lipid_class_alias_info[tmp_class]["RULE_CLASS"]
                    rgx_lst = class_rules_dct.get(rule_class, [])
    else:
        if lipid_class in class_rules_dct:
            rgx_lst = class_rules_dct.get(lipid_class, [])
        else:
            raise ValueError(f"Lipid class {lipid_class} is not supported")

    if not rgx_lst:
        rgx_lst = [r for k in class_rules_dct for r in class_rules_dct[k]]

    parsed_info_dct = get_matched_info(abbr, rgx_lst, rules_class_dct=rules_class_dct)

    if not parsed_info_dct:
        logger.warning(
            f'Can not parse abbreviation: "{abbr}", try to ignore case and try again...'
        )
        parsed_info_dct = get_matched_info(
            abbr, rgx_lst, rules_class_dct=rules_class_dct
        )
        if parsed_info_dct:
            logger.info(f'Successfully parsed abbreviation: "{abbr}"')
        else:
            logger.warning(f'Not able to parse abbreviation: "{abbr}"')
    if abbr_mod:
        parsed_info_dct["ADDITIONAL_MOD"] = abbr_mod

    return parsed_info_dct


def get_matched_info(
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
    matched_dct = {"INPUT_ABBR": abbr, "OUTPUT_INFO": {}}

    if rgx_lst:
        for rgx in rgx_lst:
            if not ignore_case:
                pass
            else:
                rgx = re.compile(rgx.pattern, flags=re.IGNORECASE)
            if rgx.match(abbr):
                rgx_match = rgx.match(abbr)
                rgx_pattern = str(rgx.pattern)
                matched_dct["OUTPUT_INFO"][rgx_pattern] = rgx_match.groupdict()

                if matched_dct["OUTPUT_INFO"][rgx_pattern]:
                    if rgx in rules_class_dct:
                        matched_dct["OUTPUT_INFO"][rgx_pattern][
                            "RULE_CLASS"
                        ] = rules_class_dct[rgx]
                    else:
                        raise ValueError(f"Can not get the lipid class of {abbr}")
                else:
                    del matched_dct["OUTPUT_INFO"][rgx_pattern]

    return matched_dct


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

    rgx_white = re.compile(r"\s*[,;\'_\-()+<>]\s*")
    abbr_rep = re.sub(rgx_white, "|", abbr)
    abbr_seg_lst = filter(None, abbr_rep.split("|"))
    for abbr_seg in abbr_seg_lst:
        seg_parsed = False
        for mod in cv_patterns_dct:
            obs_info = {}
            rgx = cv_patterns_dct[mod]
            found_segments_lst = rgx.findall(abbr_seg)
            if found_segments_lst:
                for segment_tpl in found_segments_lst:
                    segment_str = "".join(segment_tpl)
                    m = rgx.search(segment_str)
                    if m:
                        reconstruct_str = seg_to_str(m.groups(), sep="")
                        if reconstruct_str == abbr_seg:
                            obs_info = m.groupdict()
                            seg_parsed = True
            if seg_parsed and obs_info:
                cv = get_mod_cv(mod)
                obs_info["MOD"] = cv
                if cv not in mod_dct:
                    mod_dct[obs_info["MOD"]] = [obs_info]
                else:
                    mod_dct[obs_info["MOD"]].append(obs_info)
        if not seg_parsed:
            if re.match(r"O+$", abbr_seg.upper()):
                if "O" not in mod_dct:
                    mod_dct["O"] = [{"FRONT": None, "MOD": "O", "END": len(abbr_seg)}]
                else:
                    mod_dct["O"].append(
                        {"FRONT": None, "MOD": "O", "END": len(abbr_seg)}
                    )
            elif re.match(r"[aeop]", abbr_seg):
                pass
            elif len(abbr_seg) > 3:
                for mod in cv_patterns_dct:
                    rgx = cv_patterns_dct[mod]
                    found_segments_lst = rgx.findall(abbr_seg)
                    if found_segments_lst:
                        for segment_tpl in found_segments_lst:
                            segment_str = "".join(segment_tpl)
                            m = rgx.search(segment_str)
                            if m:
                                obs_info = m.groupdict()
                                cv = get_mod_cv(mod)
                                obs_info["MOD"] = cv
                                if cv not in mod_dct:
                                    mod_dct[obs_info["MOD"]] = [obs_info]
                                else:
                                    mod_dct[obs_info["MOD"]].append(obs_info)
                if "O" in mod_dct:
                    # remove replicated "O" if following modifications are detected
                    for o_mod in ["OH", "OXO", "oxo", "ep", "OO", "OOH", "Ke"]:
                        if o_mod in mod_dct:
                            mod_dct.pop("O", None)
            else:
                logger.error(
                    f"Modification {abbr} contains unsupported segment: {abbr_seg}."
                )

    return mod_dct


def remove_prefix(abbr: str = None) -> str:
    abbr = abbr.strip('"')  # remove possible side " from csv
    rgx_white = re.compile(r"\s+")
    abbr = re.sub(rgx_white, "", abbr)  # remove all white space
    if abbr.upper().startswith("OXO"):
        abbr = abbr[3:]
    if abbr.upper().startswith("OX"):
        abbr = abbr[2:]

    return abbr


def get_mod_cv(abbr: str = None) -> str:

    cv = ""
    for mod in cv_order_list:
        if mod in cv_alias_info:
            mod_alias_lst = cv_alias_info[mod]
            if abbr in mod_alias_lst:
                cv = mod

    return cv


if __name__ == "__main__":

    # s = r"O-a36:2"
    # d = parse(s)
    # print(d)

    s = r"_O_OO,OOO_2OH_oxo2_OOH"
    d = parse_mod(s)
    print(d)
