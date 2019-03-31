# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import re


class InfoOrganizer(object):

    def __init__(self):

        with open('LipidInfoStructure.json') as keys_json:
            keys_dct = json.load(keys_json)
            print(keys_dct)

    def organize_lipid(self, info_dct):
        pass


if __name__ == '__main__':

    organizer = InfoOrganizer()

