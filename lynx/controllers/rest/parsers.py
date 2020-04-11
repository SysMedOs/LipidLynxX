# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from flask_restful import reqparse

convert_get_parser = reqparse.RequestParser()
convert_get_parser.add_argument(
    "data", type=str, location=["args", "headers", "form", "json"]
)

equalizer_get_parser = reqparse.RequestParser()
equalizer_get_parser.add_argument(
    "data", type=str, location=["args", "headers", "form", "json"]
)
equalizer_get_parser.add_argument(
    "level", type=str, location=["args", "headers", "form", "json"]
)
