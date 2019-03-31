# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import LibLION.AbbrElemCalc
from LibLION.LipidNomenclature import ParserFA


def get_lpp_info(lpp_name):
    name_parser = ParserFA()
    lpp_info = name_parser.get_sum_info(lpp_name)
    print(lpp_info)


if __name__ == '__main__':

    abbr_lst = ['FA16:0', 'FA18:0', 'FA18:1', 'O-16:0', 'P-18:0']
    for abbr in abbr_lst:
        get_lpp_info(abbr)

