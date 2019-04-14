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

from LibLION.DefaultParams import logger
from LibLION.AbbrParser import AbbrParser


class Converter:

    def __init__(self, cfg: str):

        abbr_df = pd.read_excel(cfg)
        self.cfg = cfg
        self.abbr_dct = dict(zip(abbr_df['Abbreviation'].tolist(), abbr_df['epiLION'].tolist()))

    def load_file(self, file: str) -> dict:
        if os.path.isfile(file):
            if file.lower().endswith('.xlsx'):
                abbr_df = pd.read_excel(file)
            elif file.lower().endswith('.csv'):
                abbr_df = pd.read_csv(file)
            elif file.lower().endswith('.tsv'):
                abbr_df = pd.read_csv(file, sep='\t')
            else:
                abbr_df = pd.DataFrame()
                logger.error(f'Can Not load file: {file}')
        else:
            assert FileNotFoundError
            abbr_df = pd.DataFrame()
        abbr_df.fillna('', inplace=True)
        groups_lst = abbr_df.columns.tolist()
        logger.info(f'Input {len(groups_lst)} abbreviation groups: {", ".join(groups_lst)}')

        abbr_dct = {}
        if not abbr_df.empty:
            for g in groups_lst:
                col_abbr_lst = abbr_df[g].unique().tolist()
                try:
                    col_abbr_lst.remove('')
                except ValueError:
                    pass
                abbr_dct[g] = col_abbr_lst

        return abbr_dct

    def convert_table(self, input_file: str, output_file: str):

        abbr_dct = self.load_file(input_file)
        abbr_parser = AbbrParser(self.cfg)
        epilion_dct = {}
        for k in abbr_dct:
            tmp_abbr_lst = abbr_dct[k]
            epilion_lst = []
            for abbr in tmp_abbr_lst:
                abbr = re.sub(r' ', '', abbr)
                epilion_abbr = ''
                # Try to parse to software generated abbreviations
                abbr_epilion_lst = []
                if abbr_parser.is_lpptiger(abbr):
                    epilion_lpptiger_abbr = abbr_parser.parse_lpptiger(abbr)
                    logger.debug(f'{k} - LPPtiger: {abbr} -> {epilion_lpptiger_abbr}')
                    abbr_epilion_lst.append(epilion_lpptiger_abbr)
                if abbr_parser.is_lipostar(abbr)[0]:
                    epilion_lipostar_abbr = abbr_parser.parse_lipostar(abbr)
                    logger.debug(f'{k} - Lipostar: {abbr} -> {epilion_lipostar_abbr}')
                    abbr_epilion_lst.append(epilion_lipostar_abbr)
                if abbr_parser.is_lipidmaps(abbr)[0]:
                    epilion_lipidmaps_abbr = abbr_parser.parse_lipidmaps(abbr)
                    logger.debug(f'{k} - LIPIDMAPS: {abbr} -> {epilion_lipidmaps_abbr}')
                    abbr_epilion_lst.append(epilion_lipidmaps_abbr)

                if abbr_parser.is_legacy(abbr)[0]:
                    # logger.info(f'Try to parse in Legacy mode for {k} - {abbr}')
                    epilion_legacy_abbr = abbr_parser.parse_legacy(abbr)
                    logger.debug(f'{k} - Legacy: {abbr} -> {epilion_legacy_abbr}')
                    abbr_epilion_lst.append(epilion_legacy_abbr)

                if abbr_epilion_lst:
                    epilion_abbr = sorted(abbr_epilion_lst, key=len)[-1]
                else:
                    epilion_abbr = ''

                if not epilion_abbr:
                    logger.warning(f'! Can NOT convert Column {k} - {abbr}')

                epilion_lst.append(epilion_abbr)
            # logger.info(epilion_lst)
            epilion_dct[k] = epilion_lst

        # logger.info(epilion_dct)

        out_df = pd.DataFrame.from_dict(epilion_dct, orient='index').T

        if output_file.endswith('.xlsx'):
            out_df.to_excel(output_file, index=False)
        elif output_file.endswith('.csv'):
            out_df.to_csv(output_file, index=False)
        else:
            out_df.to_excel(output_file + '.xlsx', index=False)


if __name__ == '__main__':

    test_in_file = r'../Test/TestInput/test_crosscheck.xlsx'
    test_out_file = r'../Test/TestOutput/test_crosscheck_output.xlsx'
    cfg_file = r'../Configurations/LinearFA_abbreviations.xlsx'

    converter = Converter(cfg_file)

    converter.convert_table(test_in_file, test_out_file)

    logger.info('epiLion converter finished.')
