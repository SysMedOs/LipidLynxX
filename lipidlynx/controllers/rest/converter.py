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
from flask_restful import Resource
from werkzeug.utils import secure_filename

from lipidlynx.config import app_cfg_dct
from lipidlynx.config import lipidlynx_blueprint
from lipidlynx.controllers.encoder import lynx_encode
from lipidlynx.controllers.file_handler import get_table
from lipidlynx.controllers.parser import parse


@lipidlynx_blueprint.route(
    "/api/0.1/convert/name/<str:input_abbreviation>", methods=["GET", "POST"]
)
def convert_name(input_abbreviation: str):
    print("input_abbreviation", input_abbreviation)
    if input_abbreviation:
        lipidlynx_id = lynx_encode(parse(input_abbreviation))
        print("converted_ID", lipidlynx_id)
        if lipidlynx_id:
            return jsonify(
                {"code": 0, "msg": "Conversion success.", "data": lipidlynx_id}
            )
        else:
            return jsonify({"code": 1001, "msg": "Conversion failure!", "data": ""})
    else:
        return jsonify({"code": 1002, "msg": "Input error!", "data": ""})


@lipidlynx_blueprint.route(
    "/api/0.1/convert/json/<str:js_str>", methods=["GET", "POST"]
)
def convert_json(js_str: str):
    usr_dct = json.loads(js_str)
    converted_dct = {}

    if usr_dct and isinstance(usr_dct, dict):
        for k in usr_dct:
            if isinstance(k, str) and len(k) < 256:
                k_val = usr_dct[k]
                if isinstance(k_val, str):
                    in_lst = [k_val]
                elif isinstance(k_val, list):
                    in_lst = list(
                        filter(lambda x: isinstance(x, str) and len(x) > 1, k_val)
                    )
                else:
                    in_lst = []
                if in_lst:
                    output_dct = {
                        "INPUT": [],
                        "OUTPUT": [],
                        "PAIR": [],
                        "NOT_CONVERTED": [],
                    }
                    for s in in_lst:
                        lynx_id = lynx_encode(parse(s))
                        if lynx_id:
                            output_dct["INPUT"].append(s)
                            output_dct["OUTPUT"].append(lynx_id)
                            output_dct["PAIR"].append((s, lynx_id))
                        else:
                            output_dct["NOT_CONVERTED"].append(s)
                    converted_dct[k] = output_dct
    if converted_dct:
        return jsonify({"code": 0, "msg": "Json input parsed", "data": converted_dct})
    else:
        return jsonify({"code": 1003, "msg": "Json input error", "data": {}})


@lipidlynx_blueprint.route("/api/0.1/upload/<str:filename>", methods=["GET", "POST"])
def table_input(filename: str) -> str:
    abs_table_path = os.path.join(app_cfg_dct["ABS_UPLOAD_PATH"], filename)
    if os.path.isfile(abs_table_path):
        return jsonify(
            {
                "code": 0,
                "msg": "Table convert success.",
                "data": get_table(abs_table_path),
            }
        )
    else:
        return jsonify({"code": 1003, "msg": "Input error!", "data": {}})


@lipidlynx_blueprint.route("/api/0.1/download/<str:filename>", methods=["GET"])
def download(filename: str) -> str:
    if request.method == "GET":
        if os.path.isfile(os.path.join(app_cfg_dct["ABS_DOWNLOAD_PATH"], filename)):
            return send_from_directory(
                app_cfg_dct["ABS_DOWNLOAD_PATH"], filename, as_attachment=True
            )
        abort(404)
