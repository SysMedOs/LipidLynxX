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


class AbbrParser:

    def __init__(self, cfg):

        abbr_df = pd.read_excel(cfg)
        self.abbr_dct = dict(zip(abbr_df['Abbreviation'].tolist(), abbr_df['epiLION'].tolist()))

        self.fa_rgx = re.compile(r'(?P<LINK>FA|O-|P-)?(?P<NUM_C>\d{1,2})(:)(?P<NUM_DB>\d)'
                                 r'(?P<MOD_INFO>\[.*\])?(?P<REP_INFO><.*>)?')
        self.fa_lipostar_rgx = re.compile(r'(?P<LINK>FA|O-|P-)?(?P<NUM_C>\d{1,2})(:)(?P<NUM_DB>\d)'
                                          r'(?P<MOD_INFO>\(.*\))?')
        self.fa_legacy_rgx = re.compile(r'(?P<LINK1>FA|O-|P-)?(?P<NUM_C>\d{1,2})(:)(?P<NUM_DB>\d)(?P<LINK2>[ape])?'
                                        r'(?P<MOD_INFO>\([^_/]*\))?')
        self.pl_rgx = re.compile(r'(?P<LYSO>L)?(?P<PL>P[ACEGIS]|PIP[1-3]?)'
                                 r'\((?P<FA1>[^_/]*)(?P<POSITION>[_/\\])?(?P<FA2>.*)?\)(?P<HGMOD><.*>)?')
        self.pl_legacy_rgx = re.compile(r'(?P<LYSO>L)?(?P<PL1>P[ACEGIS]|PIP[1-3]?)?'
                                        r'\(?(?P<FA1>[^_/-]*)(?P<POSITION>[_/\\])?(?P<FA2>[^_/-]*)?\)?'
                                        r'(?P<PL2>-L?P[ACEGIS]|-L?PIP[1-3]?)?')

    def is_lpptiger(self, abbr: str) -> bool:

        is_correct = False

        fa_match = re.match(self.fa_rgx, abbr)
        pl_match = re.match(self.pl_rgx, abbr)

        if '[' in abbr and 'x' in abbr:
            if fa_match:
                is_correct = True
            elif pl_match:
                is_correct = True
            else:
                pass
        elif '<' in abbr:
            if fa_match:
                is_correct = True
            elif pl_match:
                is_correct = True
            else:
                pass
        else:
            pass

        return is_correct

    def parse_lpptiger(self, abbr: str) -> str:

        is_match = self.is_lpptiger(abbr)

        if is_match:
            epilion_id = abbr.replace('x', '')
        else:
            epilion_id = ''

        return epilion_id

    def is_lipostar(self, abbr: str) -> (bool, str):

        is_correct = False
        lipid_class = ''

        fa_match = re.match(self.fa_lipostar_rgx, abbr)
        pl_match = re.match(self.pl_rgx, abbr)

        if 'E' not in abbr:
            if fa_match:
                if 'a' not in abbr and 'p' not in abbr and 'o' not in abbr:
                    is_correct = True
                    lipid_class = 'FA'
            elif pl_match:
                is_correct = True
                lipid_class = 'PL'
            else:
                pass
        else:

            pass

        return is_correct, lipid_class

    def parse_lipostar_fa(self, abbr: str) -> str:
        epilion_id = ''
        fa_match = re.match(self.fa_lipostar_rgx, abbr)
        if fa_match:
            mod_lst = []
            fa_info_dct = fa_match.groupdict()
            if fa_info_dct['LINK']:
                pass
            else:
                fa_info_dct['LINK'] = 'FA'

            if int(fa_info_dct['NUM_DB']) > 0:
                mod_lst.append(f'{fa_info_dct["NUM_DB"]}DB')
            if fa_info_dct['MOD_INFO']:
                mod_lst.append(str(fa_info_dct['MOD_INFO']).strip('()'))
            if mod_lst:
                mod_str = f"[{','.join(mod_lst)}]"
            else:
                mod_str = ''

            epilion_id = f"{fa_info_dct['LINK']}{fa_info_dct['NUM_C']}:{fa_info_dct['NUM_DB']}{mod_str}"

        return epilion_id

    def parse_lipostar_pl(self, abbr: str) -> str:
        epilion_id = ''
        pl_match = re.match(self.pl_rgx, abbr)
        if pl_match:
            pl_info_dct = pl_match.groupdict()
            if pl_info_dct['LYSO']:
                pass
            else:
                pl_info_dct['LYSO'] = ''

            fa1_str = self.parse_lipostar_fa(pl_info_dct['FA1']).strip('FA')
            fa2_str = self.parse_lipostar_fa(pl_info_dct['FA2']).strip('FA')
            epilion_id = f"{pl_info_dct['LYSO']}{pl_info_dct['PL']}({fa1_str}{pl_info_dct['POSITION']}{fa2_str})"

        return epilion_id

    def parse_lipostar(self, abbr: str) -> str:

        is_match, lipid_class = self.is_lipostar(abbr)
        epilion_id = ''
        if is_match:
            if lipid_class == 'FA':
                epilion_id = self.parse_lipostar_fa(abbr)
            if lipid_class == 'PL':
                epilion_id = self.parse_lipostar_pl(abbr)

        return epilion_id

    def is_legacy(self, abbr: str) -> (bool, str):

        is_correct = False
        lipid_class = ''

        fa_match = re.match(self.fa_legacy_rgx, abbr)
        pl_match = re.match(self.pl_legacy_rgx, abbr)

        if fa_match and '/' not in abbr and '_' not in abbr:
            is_correct = True
            lipid_class = 'FA'
        else:
            if abbr in self.abbr_dct:
                is_correct = True
                lipid_class = 'FA'
            else:
                if pl_match:
                    is_correct = True
                    lipid_class = 'PL'

        return is_correct, lipid_class

    def parse_legacy_fa(self, abbr: str) -> str:
        epilion_id = ''
        fa_match = re.match(self.fa_legacy_rgx, abbr)
        if fa_match:
            mod_lst = []
            fa_info_dct = fa_match.groupdict()
            if fa_info_dct['LINK1']:
                fa_info_dct['LINK'] = fa_info_dct['LINK1']
            elif fa_info_dct['LINK1'] is None and fa_info_dct['LINK2'] == 'a':
                fa_info_dct['LINK'] = 'FA'
            elif fa_info_dct['LINK1'] is None and fa_info_dct['LINK2'] == 'p':
                fa_info_dct['LINK'] = 'P-'
            else:
                fa_info_dct['LINK'] = 'FA'

            if int(fa_info_dct['NUM_DB']) > 0:
                mod_lst.append(f'{fa_info_dct["NUM_DB"]}DB')
            if fa_info_dct['MOD_INFO']:
                mod_lst.append(str(fa_info_dct['MOD_INFO']).strip('()'))
            if mod_lst:
                mod_str = f"[{','.join(mod_lst)}]"
            else:
                mod_str = ''

            epilion_id = f"{fa_info_dct['LINK']}{fa_info_dct['NUM_C']}:{fa_info_dct['NUM_DB']}{mod_str}"
        else:
            if abbr in self.abbr_dct:
                epilion_id = self.abbr_dct[abbr]

        return epilion_id

    def parse_legacy_pl(self, abbr: str) -> str:
        epilion_id = ''
        pl_match = re.match(self.pl_legacy_rgx, abbr)
        if pl_match:
            pl_info_dct = pl_match.groupdict()
            if pl_info_dct['LYSO']:
                pass
            else:
                pl_info_dct['LYSO'] = ''

            if pl_info_dct['PL1']:
                pl_info_dct['PL'] = pl_info_dct['PL1']
            elif pl_info_dct['PL1'] is None and pl_info_dct['PL2'] is not None:
                pl_info_dct['PL'] = pl_info_dct['PL2'].strip('-')

            logger.debug(pl_info_dct)

            fa1_str = self.parse_legacy_fa(pl_info_dct['FA1'].strip('()')).strip('FA')
            fa2_str = self.parse_legacy_fa(pl_info_dct['FA2'].strip('()')).strip('FA')
            logger.debug(f"{fa2_str} {pl_info_dct['FA2']}")
            epilion_id = f"{pl_info_dct['LYSO']}{pl_info_dct['PL']}({fa1_str}{pl_info_dct['POSITION']}{fa2_str})"

        return epilion_id

    def parse_legacy(self, abbr: str) -> str:

        is_match, lipid_class = self.is_legacy(abbr)
        logger.debug(f'{is_match} {lipid_class} {abbr}')
        epilion_id = ''
        if is_match:
            if lipid_class == 'FA':
                epilion_id = self.parse_legacy_fa(abbr)
            if lipid_class == 'PL':
                epilion_id = self.parse_legacy_pl(abbr)

        return epilion_id


if __name__ == '__main__':
    pass
