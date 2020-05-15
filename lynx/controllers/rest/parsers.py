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

from flask_restful import reqparse

convert_get_parser = reqparse.RequestParser()
convert_get_parser.add_argument(
    "data", type=str, location=["args", "headers", "form", "json"]
)
convert_get_parser.add_argument(
    "export_style", type=str, location=["args", "headers", "form", "json"]
)

equalizer_get_parser = reqparse.RequestParser()
equalizer_get_parser.add_argument(
    "data", type=str, location=["args", "headers", "form", "json"]
)
equalizer_get_parser.add_argument(
    "level", type=str, location=["args", "headers", "form", "json"]
)
