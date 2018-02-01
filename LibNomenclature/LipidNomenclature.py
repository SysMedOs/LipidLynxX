# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re


class NameParserFA:

    def __init__(self):
        self.elem_dct = {'H': [1.0078250321, 0.999885],
                         'D': [2.0141017780, 0.0001157],
                         'C': [12.0, 0.9893],
                         'N': [14.0030740052, 0.99632],
                         'O': [15.9949146221, 0.99757],
                         'Na': [22.98976967, 1.0],
                         'P': [30.97376151, 1.0],
                         'S': [31.97207069, 0.9493],
                         'K': [38.9637069, 0.932581]}

        self.fa_rgx = re.compile(r'(FA)(\d{1,2})(:)(\d)')
        self.o_rgx = re.compile(r'(O-)(\d{1,2})(:)(\d)')
        self.p_rgx = re.compile(r'(P-)(\d{1,2})(:)(\d)')
        self.fa_rgx_lst = [self.fa_rgx, self.o_rgx, self.p_rgx]

    def calc_fa_mass(self, fa_info_dct):

        exactmass = (fa_info_dct['C'] * self.elem_dct['C'][0] + fa_info_dct['H'] * self.elem_dct['H'][0] +
                     fa_info_dct['O'] * self.elem_dct['O'][0])
        return round(exactmass, 6)

    def calc_fa_all_mz(self, fa_info_dct):

        exactmass = self.calc_fa_mass(fa_info_dct)
        fa_info_dct['EXACTMASS'] = exactmass
        fa_info_dct['[FA-H]-_MZ'] = round(exactmass - self.elem_dct['H'][0], 6)
        fa_info_dct['[FA-H2O-H]-_MZ'] = round(exactmass - 3 * self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
        fa_info_dct['[FA-H2O]_MZ'] = round(exactmass - 2 * self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
        fa_info_dct['[FA-H2O+H]+_MZ'] = round(exactmass - self.elem_dct['H'][0] - self.elem_dct['O'][0], 6)
        fa_info_dct['[FA-H+Na]_MZ'] = round(exactmass - self.elem_dct['H'][0] + self.elem_dct['Na'][0], 6)

        fa_info_dct['[FA-H]-_ABBR'] = '[{fa}-H]-'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O-H]-_ABBR'] = '[{fa}-H2O-H]-'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O]_ABBR'] = '[{fa}-H2O]'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H2O+H]+_ABBR'] = '[{fa}-H2O+H]+'.format(fa=fa_info_dct['ABBR'])
        fa_info_dct['[FA-H+Na]_ABBR'] = '[{fa}-H+Na]'.format(fa=fa_info_dct['ABBR'])

        return fa_info_dct

    def decode_fa(self, fa_str):

        fa_info_lst = []

        for _rgx in self.fa_rgx_lst:
            if re.match(_rgx, fa_str):
                _fa_match = re.match(_rgx, fa_str)
                fa_info_lst = _fa_match.groups()
                break

        return fa_info_lst

    def get_fa_formula(self, fa_str):

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

    def get_fa_info(self, fa_str):

        fa_info_dct = self.get_fa_formula(fa_str)
        fa_info_dct = self.calc_fa_all_mz(fa_info_dct)

        return fa_info_dct


if __name__ == '__main__':

    abbr_decoder = NameParserFA()

    abbr_lst = ['FA16:0', 'FA18:0',  'FA18:1', 'O-16:0', 'P-18:0']
    for abbr in abbr_lst:
        x = abbr_decoder.get_fa_info(abbr)
        print(x)
