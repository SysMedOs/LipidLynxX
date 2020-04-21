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
import re

import pandas as pd

from ..models.defaults import cfg_info_dct
from ..models.defaults import logger
from ..liblynx.AbbrParser import AbbrParser


class Converter:
    def __init__(self, abbr_df: pd.DataFrame = None):

        if isinstance(abbr_df, pd.DataFrame):
            self.abbr_parser = AbbrParser(abbr_df=abbr_df)
        else:
            abbr_df = pd.read_excel(cfg_info_dct["abbr_cfg"])
            self.abbr_parser = AbbrParser(abbr_df=abbr_df)
        self.abbr_dct = dict(
            zip(abbr_df["Abbreviation"].tolist(), abbr_df["LipidLynx"].tolist())
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
                logger.error(f"Can Not __load__ file: {file}")
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
        abbr_lipidlynx_lst = []
        if self.abbr_parser.is_lpptiger(abbr):
            lipidlynx_lpptiger_abbr = self.abbr_parser.parse_lpptiger(abbr)
            logger.debug(f"LPPtiger: {abbr} -> {lipidlynx_lpptiger_abbr}")
            abbr_lipidlynx_lst.append(lipidlynx_lpptiger_abbr)
        if self.abbr_parser.is_lipostar(abbr)[0]:
            lipidlynx_lipostar_abbr = self.abbr_parser.parse_lipostar(abbr)
            logger.debug(f"Lipostar: {abbr} -> {lipidlynx_lipostar_abbr}")
            abbr_lipidlynx_lst.append(lipidlynx_lipostar_abbr)
        if self.abbr_parser.is_lipidmaps(abbr)[0]:
            lipidlynx_lipidmaps_abbr = self.abbr_parser.parse_lipidmaps(abbr)
            logger.debug(f"LIPIDMAPS: {abbr} -> {lipidlynx_lipidmaps_abbr}")
            abbr_lipidlynx_lst.append(lipidlynx_lipidmaps_abbr)

        if self.abbr_parser.is_legacy(abbr)[0]:
            # logger.info(f'Try to parse in Legacy mode for {k} - {abbr}')
            lipidlynx_legacy_abbr = self.abbr_parser.parse_legacy(abbr)
            logger.debug(f"Legacy: {abbr} -> {lipidlynx_legacy_abbr}")
            abbr_lipidlynx_lst.append(lipidlynx_legacy_abbr)

        if abbr_lipidlynx_lst:
            lipidlynx_abbr = sorted(abbr_lipidlynx_lst, key=len)[-1]
        else:
            lipidlynx_abbr = ""

        if not lipidlynx_abbr:
            logger.warning(f"! Can NOT convert {abbr}")

        return lipidlynx_abbr

    def convert_table(self, input_file: str, output_file: str):

        abbr_dct = self.load_file(input_file)
        lipidlynx_dct = {}
        for k in abbr_dct:
            tmp_abbr_lst = abbr_dct[k]
            lipidlynx_lst = []
            for abbr in tmp_abbr_lst:
                lipidlynx_abbr = self.convert_abbr(abbr)
                lipidlynx_lst.append(lipidlynx_abbr)
            # logger.info(lipidlynx_lst)
            lipidlynx_dct[k] = lipidlynx_lst

        # logger.info(lipidlynx_dct)

        out_df = pd.DataFrame.from_dict(lipidlynx_dct, orient="index").T

        if output_file.endswith(".xlsx"):
            out_df.to_excel(output_file, index=False)
        elif output_file.endswith(".csv"):
            out_df.to_csv(output_file, index=False)
        else:
            out_df.to_excel(output_file + ".xlsx", index=False)

    def convert_list(self, input_list: list) -> list:

        lipidlynx_lst = []
        for abbr in input_list:
            lipidlynx_abbr = self.convert_abbr(abbr)
            lipidlynx_lst.append(lipidlynx_abbr)

        return lipidlynx_lst

    def convert_text(self, input_text: str) -> (dict, list):

        usr_abbr_lst = input_text.split("\n")

        lipidlynx_dct = {}
        bad_input_lst = []

        for abbr in usr_abbr_lst:
            lipidlynx_id = self.convert_abbr(abbr)
            if lipidlynx_id:
                lipidlynx_dct[abbr] = lipidlynx_id
            else:
                bad_input_lst.append(abbr)

        return lipidlynx_dct, bad_input_lst


if __name__ == "__main__":

    test_in_file = r"../test/test_input/test_crosscheck.xlsx"
    test_out_file = r"../test/test_output/test_crosscheck_output.xlsx"
    cfg_file = r"../configurations/defined_alias.xlsx"

    converter = Converter(cfg_file)

    converter.convert_table(test_in_file, test_out_file)

    logger.info("lynx convert_lipid finished.")
