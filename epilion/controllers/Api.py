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

from flask import Flask
from flask import render_template, redirect, url_for
from flask import request, abort, jsonify, send_from_directory
from werkzeug.utils import secure_filename
from epilion.config import DevConfig
from epilion import app_cfg_dct
from epilion import epilion_blueprint
from epilion.controllers.Encoder import lion_encode
from epilion.controllers.Parser import parse


@epilion_blueprint.route("/api/convert/<input_abbreviation>", methods=["GET", "POST"])
def convert(input_abbreviation):
    print("input_abbreviation", input_abbreviation)
    if input_abbreviation:
        epilion_id = lion_encode(parse(input_abbreviation))
        print("converted_ID", epilion_id)
        if epilion_id:
            return jsonify(
                {"code": 0, "errmsg": "Conversion success.", "result": epilion_id}
            )
        else:
            return jsonify(
                {"code": 1001, "errmsg": "Conversion failure!", "result": ""}
            )
    else:
        return jsonify({"code": 1002, "errmsg": "Input error!", "result": ""})


@epilion_blueprint.route("/api/upload", methods=["POST"])
def upload():
    usr_file = request.files["user_file"]
    if usr_file.filename:
        usr_file_name = secure_filename(usr_file.filename)
        unix_time = int(time.time())
        if usr_file_name.endswith(".csv"):
            masked_file_name = f"{unix_time}.csv"
        elif usr_file_name.endswith(".xlsx"):
            masked_file_name = f"{unix_time}.xlsx"
        elif usr_file_name.endswith(".xls"):
            masked_file_name = f"{unix_time}.xls"
        else:
            masked_file_name = ""
            abort(400, "File tye not supported.")
        usr_file.save(os.path.join(app_cfg_dct["ABS_UPLOAD_PATH"], masked_file_name))
        return jsonify(
            {"code": 0, "errmsg": "Upload success.", "fileName": masked_file_name}
        )
    else:
        return jsonify({"code": 1004, "errmsg": "Upload failure!"})


@epilion_blueprint.route("/api/download/<filename>", methods=["GET"])
def download(filename):
    if request.method == "GET":
        if os.path.isfile(os.path.join(app_cfg_dct["ABS_DOWNLOAD_PATH"], filename)):
            return send_from_directory(
                app_cfg_dct["ABS_DOWNLOAD_PATH"], filename, as_attachment=True
            )
        abort(404)
