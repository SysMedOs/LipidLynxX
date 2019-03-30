# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from functools import wraps
import re

from LibLION.DefaultParams import elem_info, elem_shifts, logger


class ParserMOD:

    def __init__(self):
        self.mod_rgx = re.compile(r'(?P<NUM>\d{1,2})(?P<MOD>\w\w)(?P<SITE_INFO>{[^\[\]{}]*})?')
        self.site_rgx = re.compile(r'(?P<SITE>\d{1,2})(?P<INFO>[EZRSab])?')

    def find_all(self, abbr: str, mode: str) -> list:
        matched_lst = []

        if mode == 'mod':
            rgx = self.mod_rgx
        elif mode == 'site':
            rgx = self.site_rgx
        else:
            rgx = self.mod_rgx

        if abbr:
            mod_match_itr = re.finditer(rgx, abbr.strip(''))
            if mod_match_itr:
                for mod_match in mod_match_itr:
                    matched_lst.append(mod_match.groupdict())
        logger.debug('matched_lst')
        logger.debug(matched_lst)

        return matched_lst

    def decode_site(self, abbr: str) -> list:
        if abbr:
            site_info_lst = self.find_all(abbr, 'site')
        else:
            site_info_lst = []

        return site_info_lst

    def decode_mod(self, abbr: str) -> list:

        fin_mod_info_lst = []

        mod_info_lst = self.find_all(abbr, 'mod')

        for mod_dct in mod_info_lst:
            if mod_dct['SITE_INFO']:
                mod_dct['SITE'] = self.decode_site(mod_dct['SITE_INFO'])
                fin_mod_info_lst.append(mod_dct)

        logger.debug(fin_mod_info_lst)

        return fin_mod_info_lst


class ParserFA:

    def __init__(self):

        self.fa_rgx = re.compile(r'(?P<LINK>FA|O-|P-)?(?P<NUM_C>\d{1,2})(:)(?P<NUM_DB>\d)'
                                 r'(?P<MOD_INFO>\[.*\])?(?P<REP_INFO><.*>)?')

        self.mod_parser = ParserMOD()

    def decode_fa(self, abbr: str) -> dict:

        fa_info_dct = {}
        if abbr:
            logger.debug(f'FA_ABBR: {abbr}')
            fa_match = re.match(self.fa_rgx, abbr)
            if fa_match:
                fa_info_dct = fa_match.groupdict()
        else:
            logger.warning('FA abbreviation is None!')

        mod_info = fa_info_dct['MOD_INFO']
        if mod_info and isinstance(mod_info, str):
            mod_info = mod_info.strip('[]')
            logger.debug(mod_info)
            mod_info_dct = self.mod_parser.decode_mod(mod_info)
        else:
            mod_info_dct = {}

        fa_info_dct['MOD'] = mod_info_dct

        return fa_info_dct

    def get_smi_fa(self, abbr: str) -> str:

        smi = None

        fa_info_dct = self.decode_fa(abbr)

        if fa_info_dct['NUM_C']:
            c_chain_lst = ['C'] * int(fa_info_dct['NUM_C'])

            if fa_info_dct['LINK']:
                if fa_info_dct['LINK'] == 'O-':
                    c_chain_lst[0] = 'OC'
                elif fa_info_dct['LINK'] == 'P-':
                    c_chain_lst[0] = r'O\C='
                    c_chain_lst[1] = r'C/'
                else:
                    c_chain_lst[0] = 'OC('
                    c_chain_lst.append(')=O')
            else:
                c_chain_lst[0] = 'OC('
                c_chain_lst.append(')=O')
        else:
            c_chain_lst = []

        if c_chain_lst:
            smi = ''.join(c_chain_lst)
        return smi

    def get_formula(self, fa_str) -> dict:

        fa_info_dct = {}
        fa_info_lst = self.decode_fa(fa_str)
        if fa_info_lst[0] == 'FA':
            fa_info_dct['ABBR'] = fa_str
            fa_info_dct['LINK'] = 'FA'
            fa_info_dct['C'] = int(fa_info_lst[1])
            fa_info_dct['H'] = 2 * int(fa_info_lst[1]) - 2 * int(fa_info_lst[3])
            fa_info_dct['O'] = 2
            fa_info_dct['DB'] = int(fa_info_lst[3])
            fa_info_dct['FORMULA'] = 'C{num_c}H{num_h}O{num_o}'.format(num_c=fa_info_dct['C'], num_h=fa_info_dct['H'],
                                                                       num_o=fa_info_dct['O'])
        elif fa_info_lst[0] == 'O-':
            fa_info_dct['ABBR'] = fa_str
            fa_info_dct['LINK'] = 'O-'
            fa_info_dct['C'] = int(fa_info_lst[1])
            fa_info_dct['H'] = 2 * int(fa_info_lst[1]) + 2 - 2 * int(fa_info_lst[3])
            fa_info_dct['O'] = 1
            fa_info_dct['DB'] = int(fa_info_lst[3])
            fa_info_dct['FORMULA'] = 'C{num_c}H{num_h}O'.format(num_c=fa_info_dct['C'], num_h=fa_info_dct['H'])
        elif fa_info_lst[0] == 'P-':
            fa_info_dct['ABBR'] = fa_str
            fa_info_dct['LINK'] = 'P-'
            fa_info_dct['C'] = int(fa_info_lst[1])
            fa_info_dct['H'] = 2 * int(fa_info_lst[1]) - 2 * int(fa_info_lst[3])
            fa_info_dct['O'] = 1
            fa_info_dct['DB'] = int(fa_info_lst[3])
            fa_info_dct['FORMULA'] = 'C{num_c}H{num_h}O'.format(num_c=fa_info_dct['C'], num_h=fa_info_dct['H'])

        return fa_info_dct

    def get_exactmass(self, lipid_info_dct: dict, decimal: int = 6) -> float:

        exactmass = 0.0

        for _elem in elem_info:
            if _elem in lipid_info_dct:
                exactmass += lipid_info_dct[_elem] * elem_info[_elem][0]
            else:
                logger.error(f'Elem{_elem} was not found')

        return round(exactmass, decimal)

    def get_mz_info(self, fa_info_dct: dict) -> dict:

        exactmass = self.get_exactmass(fa_info_dct)
        fa_info_dct['EXACTMASS'] = exactmass
        fa_info_dct['[FA-H]-_MZ'] = round(exactmass - elem_info['H'][0][0], 6)
        fa_info_dct['[FA-H2O-H]-_MZ'] = round(exactmass - 3 * elem_info['H'][0][0] - elem_info['O'][0][0], 6)
        fa_info_dct['[FA-H2O]_MZ'] = round(exactmass - 2 * elem_info['H'][0][0] - elem_info['O'][0][0], 6)
        fa_info_dct['[FA-H2O+H]+_MZ'] = round(exactmass - elem_info['H'][0][0] - elem_info['O'][0][0], 6)
        fa_info_dct['[FA-H+Na]_MZ'] = round(exactmass - elem_info['H'][0][0] + elem_info['Na'][0][0], 6)

        fa_info_dct['[FA-H]-_ABBR'] = '[{fa}-H]-'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O-H]-_ABBR'] = '[{fa}-H2O-H]-'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O]_ABBR'] = '[{fa}-H2O]'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O+H]+_ABBR'] = '[{fa}-H2O+H]+'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H+Na]_ABBR'] = '[{fa}-H+Na]'.format(fa=fa_info_dct['ABBR'])

        return fa_info_dct

    def get_sum_info(self, fa_str: str) -> dict:

        fa_info_dct = self.get_formula(fa_str)
        fa_info_dct = self.get_mz_info(fa_info_dct)

        return fa_info_dct


class ParserPL(ParserFA):

    def __init__(self):
        super(ParserPL, self).__init__()
        self.pl_rgx = re.compile(r'(?P<LYSO>L)?(?P<PL>P[ACEGIS]|PIP[1-3]?)'
                                 r'\((?P<FA1>[^_/]*)(?P<POSITION>[_/\\])?(?P<FA2>.*)?\)(?P<HGMOD><.*>)?')

    def decode_pl(self, abbr: str) -> dict:

        pl_info_dct = {}

        pl_match = re.match(self.pl_rgx, abbr)
        if pl_match:
            pl_info_dct['LIPID_INFO'] = pl_match.groupdict()

        if pl_info_dct['LIPID_INFO']:
            pl_info_dct['FA1_INFO'] = self.decode_fa(pl_info_dct['LIPID_INFO']['FA1'])
            pl_info_dct['FA2_INFO'] = self.decode_fa(pl_info_dct['LIPID_INFO']['FA2'])

        return pl_info_dct

    def get_formula_pl(self, fa_info_dct: dict, decimal: int = 6) -> float:
        pass


if __name__ == '__main__':

    fa_decoder = ParserFA()
    pl_decoder = ParserPL()

    fa_lst = ['FA18:0', '18:1', 'O-16:0', 'P-18:0',
              '20:4[4DB,2OH,1Ke]',
              '20:4[4DB{5,9,12,15},2OH{8,11},1Ke{14}]',
              '20:4[4DB{5Z,9E,12E,15E},2OH{8S,11R},1Ke{14}',
              '20:4[4DB{5Z,9E,11Z,14Z},1OH{8S}]', '9:0<CHO{@9C}>'
              ]

    for _abbr in fa_lst:
        fa = fa_decoder.decode_fa(_abbr)
        print(fa)
        _smi = fa_decoder.get_smi_fa(_abbr)
        print(_smi)

    # pl_lst = [r'PC(O-16:0/18:1)', r'PC(P-16:0_18:1)', r'PC(P-16:0/18:1)',
    #           r'PC(16:0/10:1[1DB{6E},1OH{5}]<COOH{@10C}>)']
    # for _abbr in pl_lst:
    #     pl = pl_decoder.decode_pl(_abbr)
    #     print(pl)
        # _smi = pl_decoder.get_smi_pl(_abbr)
        # print(_smi)
