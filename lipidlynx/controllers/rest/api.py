# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from flask import abort
from flask_restful import Resource, fields, marshal_with

from .parsers import convert_get_parser, convert_post_parser

from lipidlynx.models.lipid import Lipid


nested_tag_fields = {"id": fields.Integer(), "title": fields.String()}

post_fields = {"id": fields.Integer(), "title": fields.String()}


class ConvertApi(Resource):
    @marshal_with(post_fields)
    def get(self, lynx_id=None):
        if lynx_id:
            lipid = Lipid(lynx_id)

            return lipid.to_json()
