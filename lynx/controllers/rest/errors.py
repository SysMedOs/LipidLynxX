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


class ApiErrors(object):
    def __init__(self):
        # error 101x for IO related errors
        self.input_error = {"code": 1011, "msg": "Input error!", "data": ""}
        self.output_error = {"code": 1012, "msg": "Output error!", "data": ""}
        # error 102x for Converter related errors
        self.converter_error = {"code": 1020, "msg": "Converter error!", "data": ""}
        # error 103x for Converter related errors
        self.equalizer_error = {"code": 1030, "msg": "Equalizer error!", "data": ""}
