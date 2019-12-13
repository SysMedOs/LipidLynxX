# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json

from flask import abort
from flask_restful import Resource

from ..converter import convert_dict, convert_list, convert_string
from ..equalizer import Equalizer
from .errors import ApiErrors
from .parsers import convert_get_parser, equalizer_get_parser

errors = ApiErrors()


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
                return errors.input_error
        else:
            return errors.input_error


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
                return errors.output_error
        else:
            return errors.input_error


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
                return errors.output_error
        else:
            return errors.input_error


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
            return errors.input_error
        if converted_dct:
            return {"code": 0, "msg": "Json input parsed", "data": converted_dct}
        else:
            return errors.output_error


class LevelEqualizerAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/0.1/equalizer/level/
    -d 'data={"x":["PC 16:0_18:2"], "y":["PC 18:0_18:2"]}' -d 'level="D3"' -X GET
    """

    @staticmethod
    def get():
        args = equalizer_get_parser.parse_args()
        print(args)
        user_data = json.loads(args["data"])
        user_level = json.loads(args["level"])
        if isinstance(user_level, str) and user_level:
            equalizer = Equalizer(input_data=user_data, level=user_level)
            equalized_dct = equalizer.cross_match()
            if equalized_dct:
                return {"code": 0, "msg": "Json input parsed", "data": equalized_dct}
            else:
                return errors.output_error
        else:
            return errors.input_error


class MultiLevelEqualizerAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/0.1/equalizer/levels/
    -d 'data={"x":["PC 16:0_18:2"], "y":["PC 18:0_18:2"]}' -d 'level=["D3","B2"]' -X GET
    """

    @staticmethod
    def get():
        args = equalizer_get_parser.parse_args()
        print(args)
        user_data = json.loads(args["data"])
        user_level = json.loads(args["level"])
        if isinstance(user_level, list) and user_level:
            equalized_dct = {}
            for lv in user_level:
                equalizer = Equalizer(input_data=user_data, level=lv)
                equalized_dct[lv] = equalizer.cross_match()
            if equalized_dct:
                return {"code": 0, "msg": "Json input parsed", "data": equalized_dct}
            else:
                return errors.input_error
        else:
            return errors.output_error


class EqualizerAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/0.1/equalizer/
    -d 'data={"x":["PC 16:0_18:2"], "y":["PC 18:0_18:2"]}' -d 'level=Union[str, List[str]]' -X GET
    """

    @staticmethod
    def get():
        args = equalizer_get_parser.parse_args()
        print(args)
        user_data = json.loads(args["data"])
        user_level = json.loads(args["level"])
        if isinstance(user_level, str) and user_level:
            equalizer = Equalizer(input_data=user_data, level=user_level)
            equalized_dct = equalizer.cross_match()
        elif isinstance(user_level, list) and user_level:
            equalized_dct = {}
            for lv in user_level:
                equalizer = Equalizer(input_data=user_data, level=lv)
                equalized_dct[lv] = equalizer.cross_match()
        else:
            return errors.input_error
        if equalized_dct:
            return {"code": 0, "msg": "Json input parsed", "data": equalized_dct}
        else:
            return errors.input_error
