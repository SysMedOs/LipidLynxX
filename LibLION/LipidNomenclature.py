# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re

from rdkit import Chem

from LibLION.DefaultParams import elem_info, elem_shifts, logger, mod_cfg_df, pl_smi_info


class ParserMOD:

    def __init__(self):
        self.mod_rgx = re.compile(r'(?P<NUM>\d{1,2})?(?P<MOD>\w\w[\dA-M]?)(?P<SITE_INFO>{[^\[\]{}]*})?')
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
            else:
                mod_dct['SITE'] = None
            fin_mod_info_lst.append(mod_dct)

        logger.debug('MOD info:')
        logger.debug(fin_mod_info_lst)

        return fin_mod_info_lst


class ParserFA:

    def __init__(self):

        self.fa_rgx = re.compile(r'(?P<LINK>FA|O-|P-)?(?P<NUM_C>\d{1,2})(:)(?P<NUM_DB>\d)'
                                 r'(?P<MOD_INFO>\[.*\])?(?P<REP_INFO><.*>)?')

        self.mod_parser = ParserMOD()

    def is_fa(self, abbr):
        is_fa = False
        fa_match = re.match(self.fa_rgx, abbr)
        if fa_match:
            is_fa = True

        return is_fa

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
            c_chain_lst = [''] + ['C'] * int(fa_info_dct['NUM_C'])
            c_term_lst = []

            if fa_info_dct['LINK']:
                if fa_info_dct['LINK'] == 'O-':
                    c_chain_lst[0] = 'O'
                elif fa_info_dct['LINK'] == 'P-':
                    c_chain_lst[0] = r'O'
                    c_chain_lst[1] = r'\C='
                    c_chain_lst[2] = r'C/'
                else:
                    c_chain_lst[0] = 'O'
                    c_chain_lst[1] = 'C('
                    c_term_lst.append(')=O')
            else:
                c_chain_lst[0] = 'O'
                c_chain_lst[1] = 'C('
                c_term_lst.append(')=O')

            c_idx_lst = list(range(3, int(fa_info_dct['NUM_C']) + 1))
            idx_pos = 3

            for _mod in fa_info_dct['MOD']:
                _mod_code = _mod['MOD']
                if _mod['SITE']:
                    for _site in _mod['SITE']:
                        _idx = int(_site['SITE'])
                        site_code = mod_cfg_df.loc[_mod_code, 'SMI_SITE']
                        site_post_code = mod_cfg_df.loc[_mod_code, 'SMI_POST']
                        site_term_code = mod_cfg_df.loc[_mod_code, 'SMI_TERMINAL']
                        if isinstance(site_code, str):
                            c_chain_lst[_idx] = site_code
                        if isinstance(site_post_code, str):
                            c_chain_lst[_idx + 1] = site_post_code
                        if isinstance(site_term_code, str):
                            c_term_lst.append(site_term_code)
                else:
                    if _mod_code in ['DB', 'Ep']:
                        site_code = mod_cfg_df.loc[_mod_code, 'SMI_SITE']
                        site_post_code = mod_cfg_df.loc[_mod_code, 'SMI_POST']
                        _mod_count = int(_mod['NUM'])
                        _used_idx_lst = []
                        c_mod_idx = c_idx_lst[0]
                        _counter = 1
                        c_shift = 3
                        if c_shift * _mod_count > int(fa_info_dct['NUM_C']) - 2:
                            c_shift = 2  # if more C=C in chain and no bis-allylic position
                            logger.info('Too many C=C, try to remove bis-allylic positions')
                        while _counter <= _mod_count:
                            if c_mod_idx in c_idx_lst and c_mod_idx + 1 in c_idx_lst:
                                logger.info(c_mod_idx)
                                c_chain_lst[c_mod_idx] = site_code
                                c_chain_lst[c_mod_idx + 1] = site_post_code
                                _used_idx_lst.extend([c_mod_idx, c_mod_idx + 1])
                                c_idx_lst = [x for x in c_idx_lst if x not in _used_idx_lst]
                                c_mod_idx += c_shift
                                _counter += 1
                            else:
                                c_mod_idx += 1
                    else:
                        site_code = mod_cfg_df.loc[_mod_code, 'SMI_SITE']
                        site_term_code = mod_cfg_df.loc[_mod_code, 'SMI_TERMINAL']
                        _mod_count = int(_mod['NUM'])
                        _used_idx_lst = []
                        c_mod_idx = c_idx_lst[0]
                        _counter = 1
                        while _counter <= _mod_count:
                            if c_mod_idx in c_idx_lst:
                                c_chain_lst[c_mod_idx] = site_code
                                if isinstance(site_term_code, str):
                                    c_term_lst.append(site_term_code)
                                c_idx_lst.remove(c_mod_idx)
                                _counter += 1
                            else:
                                c_mod_idx += 1

        else:
            c_chain_lst = []
            c_term_lst = []

        logger.debug(c_term_lst)
        if c_term_lst:
            c_chain_lst.extend(sorted(c_term_lst, reverse=True))

        if c_chain_lst:
            smi = ''.join(c_chain_lst)

        smi = re.sub(r'\\/', r'\\', smi)
        smi = re.sub(r'/\\', r'/', smi)

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

    def is_pl(self, abbr):
        is_pl = False
        pl_match = re.match(self.pl_rgx, abbr)
        if pl_match:
            is_pl = True

        return is_pl

    def decode_pl(self, abbr: str) -> dict:

        pl_info_dct = {}

        pl_match = re.match(self.pl_rgx, abbr)
        if pl_match:
            pl_info_dct['LIPID_INFO'] = pl_match.groupdict()

        if pl_info_dct['LIPID_INFO']:
            pl_info_dct['FA1_INFO'] = self.decode_fa(pl_info_dct['LIPID_INFO']['FA1'])
            pl_info_dct['FA2_INFO'] = self.decode_fa(pl_info_dct['LIPID_INFO']['FA2'])

        pl_info_dct['LIPID_INFO']['CLASS'] = pl_info_dct['LIPID_INFO']['PL']
        if pl_info_dct['LIPID_INFO']['LYSO']:
            if isinstance(pl_info_dct['LIPID_INFO']['LYSO'], str):
                pl_info_dct['LIPID_INFO']['CLASS'] = (pl_info_dct['LIPID_INFO']['LYSO']
                                                      + pl_info_dct['LIPID_INFO']['CLASS'])

        return pl_info_dct

    def get_smi_pl(self, abbr: str) -> str:

        pl_smi = ''

        pl_info_dct = self.decode_pl(abbr)
        pl_info_dct['FA1_SMILES'] = self.get_smi_fa(pl_info_dct['LIPID_INFO']['FA1'])
        pl_info_dct['FA2_SMILES'] = self.get_smi_fa(pl_info_dct['LIPID_INFO']['FA2'])
        pl_class = pl_info_dct['LIPID_INFO']['CLASS']
        if pl_class in pl_smi_info:
            pl_hg_smi = pl_smi_info[pl_class]
            gly_part = r')C'
            pl_end = r')=O'
            # the following order is important!!
            pl_smi = ''.join([pl_hg_smi, pl_info_dct['FA2_SMILES'],
                              gly_part, pl_info_dct['FA1_SMILES'], pl_end])
        else:
            pl_smi = ''

        return pl_smi


if __name__ == '__main__':

    fa_decoder = ParserFA()
    pl_decoder = ParserPL()

    fa_lst = [
        'FA18:0', '18:1', 'O-16:0', 'P-18:0',
        '20:4[4DB,2OH,1Ke]',
        '20:4[4DB{5,9,12,15},2OH{8,11},1Ke{14}]',
        '20:4[4DB{5Z,9E,12E,15E},2OH{8S,11R},1Ke{14}]',
        '20:4[4DB{5Z,9E,11Z,14Z},1OH{8S}]', '9:0<CHO{@9C}>',
        # '20:1[PGA{8a,12b},1DB{13Z},1OH{15S}]'
    ]

    for _abbr in fa_lst:
        fa = fa_decoder.decode_fa(_abbr)
        logger.info(fa)
        _smi = fa_decoder.get_smi_fa(_abbr)
        logger.info(_abbr + ': ' + _smi)

    pl_lst = [
        r'PC(O-16:0/18:1)', r'PC(P-16:0_18:1)', r'PC(P-16:0/18:1)',
        'PC(16:0/20:4[4DB,2OH,1Ke])',
        'PC(16:0/20:4[4DB{5,9,12,15},2OH{8,11},1Ke{14}])',
        'PC(16:0/20:4[4DB{5Z,9E,12E,15E},2OH{8S,11R},1Ke{14}])',
    ]

    for _abbr in pl_lst:
        logger.info(_abbr)
        pl = pl_decoder.decode_pl(_abbr)
        logger.info(pl)
        _smi = pl_decoder.get_smi_pl(_abbr)
        logger.info(_abbr)
        logger.info(_smi)
