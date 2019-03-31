# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re


class LipidEncoder(object):

    def __init__(self):
        self.lipid_class_dct = {}
        print('Lipid Encoder initialized ...')
        self.lipid_class_dct['FA'] = {}
        self.lipid_class_dct['AL'] = {}
        self.lipid_class_dct['TG'] = {}
        self.lipid_class_dct['DG'] = {}
        self.lipid_class_dct['MG'] = {}
        self.lipid_class_dct['PA'] = {}
        self.lipid_class_dct['PC'] = {}
        self.lipid_class_dct['PE'] = {}
        self.lipid_class_dct['PG'] = {}
        self.lipid_class_dct['PI'] = {}
        self.lipid_class_dct['PS'] = {}
        self.lipid_class_dct['SM'] = {}
        self.lipid_class_dct['CR'] = {}
        self.lipid_class_dct['CL'] = {}
        self.lipid_class_dct['ST'] = {}

        self.id_level_dct = {'B': 'Bulk residues level', 'D': 'Discrete residues level',
                             'E': 'Exact position of residues level',
                             'P': 'Position specific modifications level',
                             'S': 'Stereospecific modifications level'}

    def get_class_set(self, lipid_info):
        pass


if __name__ == '__main__':
    pass
