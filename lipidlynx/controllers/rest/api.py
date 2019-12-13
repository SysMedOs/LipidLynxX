# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os
import time

from flask import request, abort, jsonify, send_from_directory
from flask_restful import Api, Resource, reqparse
from werkzeug.utils import secure_filename

from ...config import app_cfg_dct
from ...config import blueprint
from ..encoder import lynx_encode
from ..file_handler import get_table
from ..converter import convert_dict, convert_list, convert_string
from ..parser import parse
from .parsers import convert_get_parser


class StringConverterAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/0.1/converter/string/ -d 'data="PC 18:0_18:2"' -X GET
    """

    @staticmethod
    def get():
        args = convert_get_parser.parse_args()
        abbreviation = json.loads(args["data"])
        if isinstance(abbreviation, str) and abbreviation:
            converted_dct = convert_string(abbreviation)
            if converted_dct:
                return {"code": 0, "msg": "Conversion success.", "data": converted_dct}
            else:
                return {"code": 1001, "msg": "Conversion failure!", "data": ""}
        else:
            return {"code": 1002, "msg": "Input error!", "data": ""}


class ListConverterAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/0.1/converter/dict/
    -d 'data=["PC 16:0_18:2", "PC 18:0_18:2"]' -X GET
    """

    @staticmethod
    def get():
        args = convert_get_parser.parse_args()
        in_lst = json.loads(args["data"])
        if isinstance(in_lst, list) and in_lst:
            converted_dct = convert_list(in_lst)
            if converted_dct:
                return {"code": 0, "msg": "Json input parsed", "data": converted_dct}
            else:
                return {"code": 1003, "msg": "No output", "data": {}}
        else:
            return {"code": 1003, "msg": "Json input error", "data": {}}


class DictConverterAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/0.1/converter/dict/
    -d 'data={"x":["PC 16:0_18:2"], "y":["PC 18:0_18:2"]}' -X GET
    """

    @staticmethod
    def get():
        args = convert_get_parser.parse_args()
        usr_dct = json.loads(args["data"])
        if isinstance(usr_dct, dict) and usr_dct:
            converted_dct = convert_dict(usr_dct)
            if converted_dct:
                return {"code": 0, "msg": "Json input parsed", "data": converted_dct}
            else:
                return {"code": 1003, "msg": "Json input error", "data": {}}


class ConverterAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/0.1/converter/dict/
    -d 'data=Union[str|list|dict]' -X GET
    """

    @staticmethod
    def get():
        args = convert_get_parser.parse_args()
        user_data = json.loads(args["data"])
        converted_dct = {}
        if isinstance(user_data, str) and user_data:
            converted_dct = convert_string(user_data)
        elif isinstance(user_data, list) and user_data:
            converted_dct = convert_list(user_data)
        elif user_data and isinstance(user_data, dict):
            converted_dct = convert_dict(user_data)
        else:
            return {"code": 1003, "msg": "Json input error", "data": {}}
        if converted_dct:
            return {"code": 0, "msg": "Json input parsed", "data": converted_dct}
        else:
            return {"code": 1003, "msg": "Json input error", "data": {}}
