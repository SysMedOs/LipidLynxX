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

from lynx.controllers.converter import Converter
from lynx.controllers.equalizer import Equalizer
from lynx.controllers.rest.errors import ApiErrors
from lynx.controllers.rest.parsers import convert_get_parser, equalizer_get_parser

errors = ApiErrors()
converter = Converter()


def get_equalizer_params():
    args = equalizer_get_parser.parse_args()
    usr_data = json.loads(args["data"])
    if "," in args["level"]:
        try:
            usr_level = json.loads(args["level"])
        except json.decoder.JSONDecodeError:
            usr_level = args["level"].split(",")
    else:
        usr_level = args["level"]
    print(usr_level)
    return usr_data, usr_level


class StringConverterAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/{api_version}/converter/string/ -d 'data="PC 18:0_18:2"' -X GET
    or
    r=requests.get('http://127.0.0.1:5000/lipidlynx/api/{api_version}/converter/string', params='data="PC 18:0_18:2"').json()
    """

    @staticmethod
    def get():
        args = convert_get_parser.parse_args()
        abbreviation = json.loads(args["data"])
        if isinstance(abbreviation, str) and abbreviation:
            converted_dct = converter.convert_string(abbreviation)
            if converted_dct:
                return {"code": 0, "msg": "Conversion success.", "data": converted_dct}
            else:
                return errors.input_error
        else:
            return errors.input_error


class ListConverterAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/{api_version}/converter/dict/
    -d 'data=["PC 16:0_18:2", "PC 18:0_18:2"]' -X GET
    """

    @staticmethod
    def get():
        args = convert_get_parser.parse_args()
        in_lst = json.loads(args["data"])
        if isinstance(in_lst, list) and in_lst:
            converted_dct = converter.convert_list(in_lst)
            if converted_dct:
                return {"code": 0, "msg": "Json input parsed", "data": converted_dct}
            else:
                return errors.output_error
        else:
            return errors.input_error


class DictConverterAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/{api_version}/converter/dict/
    -d 'data={"x":["PC 16:0_18:2"], "y":["PC 18:0_18:2"]}' -X GET
    """

    @staticmethod
    def get():
        args = convert_get_parser.parse_args()
        usr_dct = json.loads(args["data"])
        if isinstance(usr_dct, dict) and usr_dct:
            converted_dct = converter.convert_dict(usr_dct)
            if converted_dct:
                return {"code": 0, "msg": "Json input parsed", "data": converted_dct}
            else:
                return errors.output_error
        else:
            return errors.input_error


class ConverterAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/{api_version}/converter/dict/
    -d 'data=Union[str|list|dict]' -X GET
    """

    @staticmethod
    def get():
        args = convert_get_parser.parse_args()
        usr_data = json.loads(args["data"])
        converted_dct = {}
        if isinstance(usr_data, str) and usr_data:
            converted_dct = converter.convert_string(usr_data)
        elif isinstance(usr_data, list) and usr_data:
            converted_dct = converter.convert_list(usr_data)
        elif usr_data and isinstance(usr_data, dict):
            converted_dct = converter.convert_dict(usr_data)
        else:
            return errors.input_error
        if converted_dct:
            return {"code": 0, "msg": "Json input parsed", "data": converted_dct}
        else:
            return errors.output_error


class LevelEqualizerAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/{api_version}/equalizer/level/ -d 'data={"x":["PC 16:0_18:2"], "y":["PC 18:0_18:2"]}' -d 'level="D3"' -X GET

    r=requests.get('http://127.0.0.1:5000/lipidlynx/api/{api_version}/equalizer',
    params={"data":'{"x":["PC 16:0_18:2"],"y":["PC 18:0_18:2"]}', "level":"D3"}).json()
    """

    @staticmethod
    def get():
        usr_data, usr_level = get_equalizer_params()
        if isinstance(usr_level, str) and usr_level:
            equalizer = Equalizer(input_data=usr_data, level=usr_level)
            equalized_dct = equalizer.cross_match()
            if equalized_dct:
                return {"code": 0, "msg": "Json input parsed", "data": equalized_dct}
            else:
                return errors.output_error
        else:
            return errors.input_error


class MultiLevelEqualizerAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/{api_version}/equalizer/levels/
    -d 'data={"x":["PC 16:0_18:2"], "y":["PC 18:0_18:2"]}' -d 'level=["D3","B2"]' -X GET

    r=requests.get('http://127.0.0.1:5000/lipidlynx/api/{api_version}/equalizer',
    params={"data":'{"x":["PC 16:0_18:2"],"y":["PC 18:0_18:2"]}', "level":'["D3","B3"]'}).json()
    """

    @staticmethod
    def get():
        usr_data, usr_level = get_equalizer_params()
        if isinstance(usr_level, list) and usr_level:
            equalized_dct = {}
            for lv in usr_level:
                equalizer = Equalizer(input_data=usr_data, level=lv)
                equalized_dct[lv] = equalizer.cross_match()
            if equalized_dct:
                return {"code": 0, "msg": "Json input parsed", "data": equalized_dct}
            else:
                return errors.input_error
        else:
            return errors.output_error


class EqualizerAPI(Resource):
    """
    $ curl http://127.0.0.1:5000/lipidlynx/api/{api_version}/equalizer/
    -d 'data={"x":["PC 16:0_18:2"], "y":["PC 18:0_18:2"]}' -d 'level=Union[str, List[str]]' -X GET

    or

    r=requests.get('http://127.0.0.1:5000/lipidlynx/api/{api_version}/equalizer',
    params={"data":'{"x":["PC 16:0_18:2"],"y":["PC 18:0_18:2"]}', "level":'["D3","B3"]'}).json()
    """

    @staticmethod
    def get():
        usr_data, usr_level = get_equalizer_params()
        equalized_dct = {}
        if isinstance(usr_level, str) and usr_level:
            if "[" in usr_level:
                usr_level = json.loads(usr_level)
            else:
                usr_level = [usr_level]
        elif isinstance(usr_level, list) and usr_level:
            pass
        else:
            return errors.input_error
        for lv in usr_level:
            equalizer = Equalizer(input_data=usr_data, level=lv)
            equalized_dct[lv] = equalizer.cross_match()

        if equalized_dct:
            return {"code": 0, "msg": "Json input parsed", "data": equalized_dct}
        else:
            return errors.input_error
