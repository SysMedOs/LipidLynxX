# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List, Tuple

import pandas as pd

from epilion.libLION.DefaultParams import logger


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

    class_rgx_dct = {}
    rgx_class_dct = {}

    for i, r in rules_df.iterrows():
        if isinstance(r["CLASS"], str):
            rgx_str = r["REGULAR_EXPRESSION"].strip('"')
            rgx = re.compile(rgx_str)

            rgx_checker = True
            if isinstance(r["EXAMPLE"], str):
                if not rgx.match(r["EXAMPLE"]):
                    rgx_checker = False
                    logger.warning(
                        f'Rule example: "{r["EXAMPLE"]}" NOT fit with rule: "{rgx_str}" -> skipped...'
                    )

            if rgx_checker:
                rgx_class_dct[rgx] = r["CLASS"]
                if r["CLASS"] not in class_rgx_dct:
                    class_rgx_dct[r["CLASS"]] = [rgx]
                else:
                    class_rgx_dct[r["CLASS"]].append(rgx)

    return class_rgx_dct, rgx_class_dct


def parse(abbr: str, class_rgx_dct: dict, rgx_class_dct: dict, lipid_class: str = None) -> dict:

    parsed_info_dct = {}
    rgx_lst = []
    detected_lipid_class = None

    if not lipid_class:
        if abbr.upper().startswith(('FA', 'O-', 'P-')):
            rgx_lst = class_rgx_dct.get('FA', [])
        elif abbr.upper().startswith(('PL', 'PA', 'PC', 'PE', 'PG', 'PS', 'PI')):
            rgx_lst = class_rgx_dct.get('PL', [])
        elif abbr.upper().startswith(('SM', 'SP')):
            rgx_lst = class_rgx_dct.get('SM', [])
        elif abbr.upper().startswith(('CER')):
            rgx_lst = class_rgx_dct.get('Cer', [])
        elif abbr.upper().startswith(('MG', 'MAG')):
            rgx_lst = class_rgx_dct.get('MG', [])
        elif abbr.upper().startswith(('DG', 'DAG')):
            rgx_lst = class_rgx_dct.get('DG', [])
        elif abbr.upper().startswith(('TG', 'TAG')):
            rgx_lst = class_rgx_dct.get('TG', [])
    else:
        if lipid_class in class_rgx_dct:
            rgx_lst = class_rgx_dct.get(lipid_class, [])

    if not rgx_lst:
        rgx_lst = [r for k in class_rgx_dct for r in class_rgx_dct[k]]

    print(abbr, rgx_lst, len(rgx_lst))
    if rgx_lst:
        for rgx in rgx_lst:
            if rgx.match(abbr):
                rgx_match = rgx.match(abbr)
                parsed_info_dct = rgx_match.groupdict()
                if parsed_info_dct and parsed_info_dct.get('CLASS', None) is None:
                    if rgx in rgx_class_dct:
                        parsed_info_dct['CLASS'] = rgx_class_dct[rgx]

    if not parsed_info_dct:
        logger.warning(f'Can not parse abbreviation: "{abbr}", try to ignore case and try again...')
        for rgx in rgx_lst:
            rgx_no_case = re.compile(rgx.pattern, flags=re.IGNORECASE)
            if rgx_no_case.match(abbr):
                rgx_match = rgx_no_case.match(abbr)
                parsed_info_dct = rgx_match.groupdict()
                if parsed_info_dct and parsed_info_dct.get('CLASS', None) is None:
                    if rgx in rgx_class_dct:
                        parsed_info_dct['CLASS'] = rgx_class_dct[rgx]
        if parsed_info_dct:
            logger.info(f'Successfully parsed abbreviation: "{abbr}"')
        else:
            logger.warning(f'Notable to parse abbreviation: "{abbr}"')

    return parsed_info_dct


if __name__ == "__main__":
    usr_rules_file = r"../configurations/rules.csv"

    usr_class_rgx_dct, usr_rgx_class_dct = build_parser(usr_rules_file)
    print(usr_class_rgx_dct)
    usr_abbr_lst = [r'FA 20:4;O2', r'fa 20:4;O2', r'Test']

    for usr_abbr in usr_abbr_lst:
        usr_parsed_info_dct = parse(usr_abbr, usr_class_rgx_dct, usr_rgx_class_dct)
        print(usr_parsed_info_dct)
