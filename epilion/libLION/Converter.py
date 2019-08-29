# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os
import re

import pandas as pd

from epilion.libLION.DefaultParams import logger
from epilion.libLION.AbbrParser import AbbrParser


class Converter:
    def __init__(self, cfg: str = None, abbr_df: pd.DataFrame = None):

        if isinstance(abbr_df, pd.DataFrame):
            self.abbr_parser = AbbrParser(abbr_df=abbr_df)
        else:
            abbr_df = pd.read_excel(cfg)
            self.abbr_parser = AbbrParser(cfg=cfg)
        self.abbr_dct = dict(
            zip(abbr_df["Abbreviation"].tolist(), abbr_df["epiLION"].tolist())
        )

    @staticmethod
    def load_file(file: str) -> dict:

        if os.path.isfile(file):
            if file.lower().endswith(".xlsx"):
                abbr_df = pd.read_excel(file)
            elif file.lower().endswith(".csv"):
                abbr_df = pd.read_csv(file)
            elif file.lower().endswith(".tsv"):
                abbr_df = pd.read_csv(file, sep="\t")
            else:
                abbr_df = pd.DataFrame()
                logger.error(f"Can Not load file: {file}")
        else:
            raise FileNotFoundError
        abbr_df.fillna("", inplace=True)
        groups_lst = abbr_df.columns.tolist()
        logger.info(
            f'Input {len(groups_lst)} abbreviation groups: {", ".join(groups_lst)}'
        )

        abbr_dct = {}
        if not abbr_df.empty:
            for g in groups_lst:
                col_abbr_lst = abbr_df[g].unique().tolist()
                try:
                    col_abbr_lst.remove("")
                except ValueError:
                    pass
                abbr_dct[g] = col_abbr_lst

        return abbr_dct

    def convert_abbr(self, abbr: str) -> str:

        abbr = re.sub(r" ", "", abbr)
        # Try to parse to software generated abbreviations
        abbr_epilion_lst = []
        if self.abbr_parser.is_lpptiger(abbr):
            epilion_lpptiger_abbr = self.abbr_parser.parse_lpptiger(abbr)
            logger.debug(f"LPPtiger: {abbr} -> {epilion_lpptiger_abbr}")
            abbr_epilion_lst.append(epilion_lpptiger_abbr)
        if self.abbr_parser.is_lipostar(abbr)[0]:
            epilion_lipostar_abbr = self.abbr_parser.parse_lipostar(abbr)
            logger.debug(f"Lipostar: {abbr} -> {epilion_lipostar_abbr}")
            abbr_epilion_lst.append(epilion_lipostar_abbr)
        if self.abbr_parser.is_lipidmaps(abbr)[0]:
            epilion_lipidmaps_abbr = self.abbr_parser.parse_lipidmaps(abbr)
            logger.debug(f"LIPIDMAPS: {abbr} -> {epilion_lipidmaps_abbr}")
            abbr_epilion_lst.append(epilion_lipidmaps_abbr)

        if self.abbr_parser.is_legacy(abbr)[0]:
            # logger.info(f'Try to parse in Legacy mode for {k} - {abbr}')
            epilion_legacy_abbr = self.abbr_parser.parse_legacy(abbr)
            logger.debug(f"Legacy: {abbr} -> {epilion_legacy_abbr}")
            abbr_epilion_lst.append(epilion_legacy_abbr)

        if abbr_epilion_lst:
            epilion_abbr = sorted(abbr_epilion_lst, key=len)[-1]
        else:
            epilion_abbr = ""

        if not epilion_abbr:
            logger.warning(f"! Can NOT convert {abbr}")

        return epilion_abbr

    def convert_table(self, input_file: str, output_file: str):

        abbr_dct = self.load_file(input_file)
        epilion_dct = {}
        for k in abbr_dct:
            tmp_abbr_lst = abbr_dct[k]
            epilion_lst = []
            for abbr in tmp_abbr_lst:
                epilion_abbr = self.convert_abbr(abbr)
                epilion_lst.append(epilion_abbr)
            # logger.info(epilion_lst)
            epilion_dct[k] = epilion_lst

        # logger.info(epilion_dct)

        out_df = pd.DataFrame.from_dict(epilion_dct, orient="index").T

        if output_file.endswith(".xlsx"):
            out_df.to_excel(output_file, index=False)
        elif output_file.endswith(".csv"):
            out_df.to_csv(output_file, index=False)
        else:
            out_df.to_excel(output_file + ".xlsx", index=False)

    def convert_list(self, input_list: list) -> list:

        epilion_lst = []
        for abbr in input_list:
            epilion_abbr = self.convert_abbr(abbr)
            epilion_lst.append(epilion_abbr)

        return epilion_lst

    def convert_text(self, input_text: str) -> (dict, list):

        usr_abbr_lst = input_text.split("\n")

        epilion_dct = {}
        bad_input_lst = []

        for abbr in usr_abbr_lst:
            epilion_id = self.convert_abbr(abbr)
            if epilion_id:
                epilion_dct[abbr] = epilion_id
            else:
                bad_input_lst.append(abbr)

        return epilion_dct, bad_input_lst


if __name__ == "__main__":

    test_in_file = r"../test/TestInput/test_crosscheck.xlsx"
    test_out_file = r"../test/TestOutput/test_crosscheck_output.xlsx"
    cfg_file = r"../configurations/LinearFA_abbreviations.xlsx"

    converter = Converter(cfg_file)

    converter.convert_table(test_in_file, test_out_file)

    logger.info("epiLion converter finished.")
