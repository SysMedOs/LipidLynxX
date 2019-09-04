# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import configparser
import re
from typing import Tuple, Dict

import pandas as pd

from epilion.controllers.GeneralFunctions import get_abs_path
from epilion.controllers.Logger import logger


def load_cfg_info(cfg_path: str = None) -> Dict[str, str]:
    cfg_path_dct = {}
    default_fields = ["cv", "rules", "mod_cfg", "abbr_cfg"]
    config = configparser.ConfigParser()
    if cfg_path and isinstance(cfg_path, str):
        config_path = get_abs_path(cfg_path)
    else:
        try:
            config_path = get_abs_path("config.ini")
        except FileNotFoundError:
            config_path = get_abs_path("configure.ini")

    config.read(config_path)
    if config.has_section("settings"):
        user_cfg = "settings"
    elif config.has_section("default"):
        user_cfg = "default"
    else:
        user_cfg = ""
        raise ValueError(f"Can not load settings from file {config_path}")

    if len(user_cfg) > 2:
        options = config.options(user_cfg)
        for field in default_fields:
            if field in options:
                cfg_path_dct[field] = get_abs_path(config.get(user_cfg, field))

    return cfg_path_dct


def build_parser(rules_file: str) -> Tuple[dict, dict]:
    """
    Read predefined rules from configurations folder and export as a dictionary
    Args:
        rules_file (str): the path for the rules file

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

    return class_rules_dct, rules_class_dct


def get_cv_lst(cv_file: str) -> list:
    cv_df = pd.read_csv(cv_file)
    cv_lst = cv_df["CV"].tolist()

    return cv_lst


def build_mod_parser(cv_file: str) -> dict:
    cv_lst = get_cv_lst(cv_file)
    cv_patterns_dct = {}
    for cv in cv_lst:
        cv_patterns_dct[cv] = re.compile(
            r"(?P<FRONT>\d\d?)?(?P<MOD>{mod})(?P<END>\d\d?)?".format(mod=cv)
        )

    return cv_patterns_dct
