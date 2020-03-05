# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List, Union

from lynx.models.log import logger
from lynx.models.defaults import (
    class_rgx_dct,
    rgx_class_dct,
    cv_rgx_dct,
    lipid_class_alias_info,
    cv_order_list,
    cv_alias_info,
)
from lynx.controllers.general_functions import seg_to_str


def parse(lipid_name: str, rules: dict) -> Dict[str, Union[str, dict]]:

    """
    Main parser to read input abbreviations
    Args:
        lipid_name: input lipid abbreviation to be converted
        rules: the predefined dict in the form of lipid_class: re.compile(r"rule")

    Returns:
        parsed_info_dct: parsed information stored as dict

    """

    parsed_info_dct = {}

    for c in rules:
        c_search_rgx = rules[c].get("SEARCH", None)
        c_match_rgx_dct = rules[c].get("MATCH", None)
        if isinstance(c_search_rgx, re.Pattern) and isinstance(c_match_rgx_dct, dict):
            class_search = c_search_rgx.search(lipid_name)
            matched_info_dct = {}
            if class_search:
                for m in c_match_rgx_dct:
                    m_pattern = c_match_rgx_dct[m]["MATCH"]  # type: re.Pattern
                    m_groups = c_match_rgx_dct[m]["GROUPS"]  # type: list
                    m_match = m_pattern.match(lipid_name)
                    if m_match:
                        matched_dct = {}
                        matched_groups = m_match.groupdict()
                        for g in m_groups:
                            matched_dct[g] = matched_groups.get(g, "")
                        matched_info_dct[m] = matched_dct
                parsed_info_dct[c] = matched_info_dct
            else:
                pass
        else:
            logger.error(f"Cannot load rules correctly: {rules}")

    if not parsed_info_dct:
        logger.error(f"Failed to decode Lipid: {lipid_name}")

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
        rgx_lst: list of possible regular expression patterns for in the form of List[re.compile]
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
    examples = [
        "LPE 16:0",
        "LPE 16:0/0:0",
        "LPE O-18:1",
        "LPE P-16:0",
        "PE 34:2",
        "PE 16:0_18:2",
        "PE O-34:2",
        "PE O-16:0_18:2",
        "PE P-34:2",
        "PE P-16:0_18:2",
        "PIP 34:1",
        "PIP 16:0_18:1",
        "LPIP 16:0",
        "LPIP 16:0/0:0",
        "PIP2 34:2",
        "PIP2 16:0_18:2",
        "LPIP2 16:0",
        "LPIP2 16:0/0:0",
        "PIP3 34:2",
        "PIP3 16:0_18:2",
        "LPIP3 16:0",
        "LPIP3 16:0/0:0",
        "PEtOH 34:2",
        "PEtOH 16:0_18:2",
        "BMP 34:1",
        "BMP 16:0_18:1",
        "MG(16:0)",
        "MG(0:0/0:0/16:0)",
        "MG(0:0/16:0/0:0)",
        "MG(16:0/0:0/0:0)",
        "DG(34:2)",
        "DG(16:0_18:2)",
        "DG(O-34:2)",
        "DG(O-16:0_18:2)",
        "DG(P-34:2)",
        "DG(P-16:0_18:2)",
        "DG(P-16:0/0:0/18:2)",
        "TG(52:2)",
        "TG(16:0_18:0_18:2)",
        "TG(16:0/18:2/18:0)",
        "TG(O-52:2)",
        "TG(O-16:0_18:0_18:2)",
        "TG(O-16:0/18:2/18:0)",
        "TG(P-52:2)",
        "TG(P-16:0_18:0_18:2)",
        "TG(P-16:0/18:2/18:0)",
    ]

    js_folder = r"../configurations/rules/input"

    from lynx.controllers.rules_reader import build_all_rules

    all_rules = build_all_rules(js_folder)

    for e in examples:
        p = parse(e, rules=all_rules)
        logger.debug(e)
        logger.debug(p)

    logger.info("fin")
