# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from flask_restful import reqparse

convert_get_parser = reqparse.RequestParser()
convert_get_parser.add_argument('id', type=str, location=['args', 'headers', 'form', 'json'])

convert_post_parser = reqparse.RequestParser()
convert_post_parser.add_argument(
    'lipid_name',
    type=str,
    required=True,
    help="Input lipid name is required"
)
