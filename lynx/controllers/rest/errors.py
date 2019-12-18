# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
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
