# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# LipidLynxX is Dual-licensed
#   For academic and non-commercial use: GPLv2 License:
#   For commercial use: please contact the SysMedOs team by email.
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: lipid annotations converter for large scale lipidomics and epilipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os
from typing import Dict, List, Tuple

from natsort import natsorted
import pandas as pd
import regex as re

from lynx.models.rules import InputRules, OutputRules
from lynx.utils.file_handler import load_folder
from lynx.utils.log import app_logger


def build_parser(rules_file: str, logger=app_logger) -> Tuple[dict, dict]:
    """
    Read predefined rules from configurations folder and export as a dictionary

    Args:
        rules_file: the path for the rules file
        logger: for log output

    Returns:
        dict contains the regular expressions as a dict

    """

    rules_df = pd.read_csv(
        rules_file,
        header=0,
        index_col=False,
        skip_blank_lines=True,
        na_filter=True,
        na_values=None,
    )
    rules_df.drop_duplicates(
        subset=["CLASS", "REGULAR_EXPRESSION"], keep="first", inplace=True
    )

    class_rules_dct = {}
    rules_class_dct = {}

    for i, r in rules_df.iterrows():
        if isinstance(r["CLASS"], str):
            rules_str = r["REGULAR_EXPRESSION"].strip('"')
            rules = re.compile(rules_str)

            rules_checker = True
            if isinstance(r["EXAMPLE"], str):
                if not rules.match(r["EXAMPLE"]):
                    rules_checker = False
                    logger.warning(
                        f'Rule example: "{r["EXAMPLE"]}" NOT fit with rule: "{rules_str}" -> skipped...'
                    )

            if rules_checker:
                rules_class_dct[rules] = r["CLASS"]
                if r["CLASS"] not in class_rules_dct:
                    class_rules_dct[r["CLASS"]] = [rules]
                else:
                    class_rules_dct[r["CLASS"]].append(rules)

                logger.debug(
                    f'Rule added: "{r["CLASS"]}" -> "{r["REMARK"]}" -> "{r["EXAMPLE"]}"'
                )

    return class_rules_dct, rules_class_dct


def get_cv_lst(cv_file: str) -> list:
    cv_df = pd.read_csv(cv_file)
    cv_lst = cv_df["CV"].tolist()

    return cv_lst


def build_mod_parser(cv_alias_info: Dict[str, List[str]]) -> dict:
    cv_patterns_dct = {}
    for cv in cv_alias_info:
        alias_lst = cv_alias_info[cv]
        for alia in alias_lst:
            cv_patterns_dct[alia] = re.compile(
                r"(\s*[;_+]\s*)?(?P<FRONT>\d\d?[xX]?)?(?P<MOD>{mod})(?P<END>\d\d?)?(?P<REPLACE>@[CHNOP]\d\d?)?".format(
                    mod=alia
                )
            )

    return cv_patterns_dct


def build_input_rules(folder: str, logger=app_logger) -> dict:

    input_rules = {}
    file_path_lst = load_folder(folder, file_type=".json")
    # logger.debug(f"Fund JSON config files: \n {file_path_lst}")

    for f in file_path_lst:
        temp_rules = InputRules(f, logger)
        idx_lst = [os.path.basename(f)] + temp_rules.sources
        idx = "#".join(idx_lst)
        for c in temp_rules.rules:
            c_class_str = temp_rules.rules[c].get("CLASS", "")
            c_lmsd_classes = temp_rules.rules[c].get("LMSD_CLASSES", [])
            existed_c_info = input_rules.get(c_class_str, {})
            c_rgx = existed_c_info.get("SEARCH", re.compile(c_class_str))
            c_pattern = existed_c_info.get("MATCH", {})
            c_lmsd_classes.extend(existed_c_info.get("LMSD_CLASSES", []))
            c_lmsd_classes = natsorted(list(set(c_lmsd_classes)))
            c_info = temp_rules.rules[c]
            del c_info["CLASS"]
            del c_info["LMSD_CLASSES"]
            c_pattern[idx] = c_info

            input_rules[c_class_str] = {
                "LMSD_CLASSES": c_lmsd_classes,
                "SEARCH": c_rgx,
                "MATCH": c_pattern,
                "RESIDUES_SEPARATOR": temp_rules.separators.get(
                    "RESIDUES_SEPARATOR", "_|/"
                ),
                "SEPARATOR_LEVELS": temp_rules.separators.get("SEPARATOR_LEVELS", {}),
                "MAX_RESIDUES": temp_rules.rules[c].get("MAX_RESIDUES", 1),
            }

    # logger.debug(input_rules)

    return input_rules


def build_output_rules(folder: str, logger=app_logger) -> dict:

    output_rules = {}
    file_path_lst = load_folder(folder, file_type=".json")
    # logger.debug(f"Fund JSON config files: \n {file_path_lst}")

    for f in file_path_lst:
        temp_rules = OutputRules(f)
        idx = f"{temp_rules.nomenclature}@{temp_rules.date}"
        output_rules[idx] = temp_rules.rules

    logger.debug(output_rules)

    return output_rules


def load_output_rule(output_rules: dict, rule: str = "LipidLynxX"):
    output_rule_info = output_rules.get(rule, None)
    if not output_rule_info:
        rule_ver = 0
        for o_rule in output_rules:
            o_rule_lst = o_rule.split("@")
            if len(o_rule_lst) == 2:
                nomenclature = o_rule_lst[0]
                version = int(o_rule_lst[1])
                if rule.lower().startswith(nomenclature.lower()) and version > rule_ver:
                    output_rule_info = output_rules.get(o_rule, None)
    else:
        raise ValueError(f"Cannot load output rule: {rule}")
    if output_rule_info:
        return output_rule_info
    else:
        raise ValueError(f"Cannot load output rule: {rule}")
