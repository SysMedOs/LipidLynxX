# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List

import pandas as pd

from epilion.libLION.DefaultParams import logger


def build_parser(rules_file: str) -> Dict[str, List[re.compile]]:
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

    rgx_dct = {}

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
                if r["CLASS"] not in rgx_dct:
                    rgx_dct[r["CLASS"]] = [rgx]
                else:
                    rgx_dct[r["CLASS"]].append(rgx)

    return rgx_dct


def parse(abbr:str, rgx_dct: Dict[str, List[re.compile]]) -> dict:

    parsed_info_dct = {}

    if abbr.startswith(('FA', 'O-', 'P-')):
        rgx_lst = rgx_dct.get('FA', [])
        if rgx_lst:
            for rgx in rgx_lst:
                if rgx.match(abbr):
                    rgx_match = rgx.match(abbr)
                    parsed_info_dct = rgx_match.groupdict()

    return parsed_info_dct


if __name__ == "__main__":
    usr_rules_file = r"../configurations/rules.csv"

    usr_rgx_dct = build_parser(usr_rules_file)
    print(usr_rgx_dct)
    usr_abbr_lst = [r'FA 20:4;O2']

    for usr_abbr in usr_abbr_lst:
        usr_parsed_info_dct = parse(usr_abbr, usr_rgx_dct)
        print(usr_parsed_info_dct)
