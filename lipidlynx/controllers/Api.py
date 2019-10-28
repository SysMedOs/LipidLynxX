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
from werkzeug.utils import secure_filename

from lipidlynx import app_cfg_dct
from lipidlynx import lipidlynx_blueprint
from lipidlynx.controllers.Encoder import lynx_encode
from lipidlynx.controllers.FileIO import get_table
from lipidlynx.controllers.Parser import parse


@lipidlynx_blueprint.route(
    "/api/single_convert/<input_abbreviation>", methods=["GET", "POST"]
)
def single_convert(input_abbreviation):
    print("input_abbreviation", input_abbreviation)
    if input_abbreviation:
        lipidlynx_id = lynx_encode(parse(input_abbreviation))
        print("converted_ID", lipidlynx_id)
        if lipidlynx_id:
            return jsonify(
                {"code": 0, "errmsg": "Conversion success.", "result": lipidlynx_id}
            )
        else:
            return jsonify(
                {"code": 1001, "errmsg": "Conversion failure!", "result": ""}
            )
    else:
        return jsonify({"code": 1002, "errmsg": "Input error!", "result": ""})


@lipidlynx_blueprint.route("/api/table_input/<filename>", methods=["GET", "POST"])
def table_input(filename):
    abs_table_path = os.path.join(app_cfg_dct["ABS_UPLOAD_PATH"], filename)
    if os.path.isfile(abs_table_path):
        return jsonify(
            {
                "code": 0,
                "errmsg": "Table convert success.",
                "result": get_table(abs_table_path),
            }
        )
    else:
        return jsonify(
            {
                "code": 1003,
                "errmsg": "Input error!",
                "result": get_table(abs_table_path),
            }
        )


@lipidlynx_blueprint.route("/api/download/<filename>", methods=["GET"])
def download(filename):
    if request.method == "GET":
        if os.path.isfile(os.path.join(app_cfg_dct["ABS_DOWNLOAD_PATH"], filename)):
            return send_from_directory(
                app_cfg_dct["ABS_DOWNLOAD_PATH"], filename, as_attachment=True
            )
        abort(404)
