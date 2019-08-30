# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
import os
import time

from flask import Flask
from flask import Blueprint
from flask import render_template, redirect, url_for
from flask import request, abort, jsonify, send_from_directory
from flask_wtf import Form
from werkzeug.utils import secure_filename
from wtforms import StringField, TextAreaField
from wtforms.validators import DataRequired, Length

from epilion.config import DevConfig
from epilion.controllers.AbbrConverter import web_converter
from epilion.controllers.LionParser import parse_epilion


class ConverterInputForm(Form):

    input_id_str = TextAreaField(
        "Paste lipid abbreviations here:", validators=[DataRequired()]
    )


class ParserInputForm(Form):

    lion_id_str = StringField(
        "Paste one lipid abbreviation here:",
        validators=[DataRequired(), Length(max=255)],
    )


epilion_blueprint = Blueprint(
    "epilion",
    __name__,
    template_folder=r"epilion/templates",
    static_folder=r"epilion/static",
    url_prefix="/epilion",
)
app_cfg_dct = {
    "ABS_BASE_PATH": os.path.abspath(os.path.dirname(__file__)),
    "UPLOAD_FOLDER": os.path.join("epilion", "uploads"),
    "DOWNLOAD_FOLDER": os.path.join("epilion", "downloads"),
}

app_cfg_dct["ABS_UPLOAD_PATH"] = os.path.join(
    app_cfg_dct["ABS_BASE_PATH"], app_cfg_dct["UPLOAD_FOLDER"]
)
app_cfg_dct["ABS_DOWNLOAD_PATH"] = os.path.join(
    app_cfg_dct["ABS_BASE_PATH"], app_cfg_dct["DOWNLOAD_FOLDER"]
)

app = Flask(__name__)
app.config.from_object(DevConfig)


@app.route("/")
def index():
    return redirect(url_for("epilion.home"))


@epilion_blueprint.route("/")
def home():
    return render_template("home.html")


@epilion_blueprint.route("/converter", methods=("GET", "POST"))
def converter():
    conv_in_form = ConverterInputForm()
    submitted = 0
    if conv_in_form.validate_on_submit():
        out_dct, bad_in_lst = web_converter(conv_in_form.input_id_str.data)
        submitted = 1
    else:
        out_dct = {}
        bad_in_lst = []

    if not out_dct:
        bad_in_lst = [""]

    output_name = "test_converter_output.csv"

    return render_template(
        "converter.html",
        out_dct=out_dct,
        bad_in_lst=bad_in_lst,
        in_form=conv_in_form,
        submitted=submitted,
        output_name=output_name,
    )


@epilion_blueprint.route("/parser", methods=("GET", "POST"))
def parser():
    in_form = ParserInputForm()
    if in_form.validate_on_submit():
        out_dct = parse_epilion(in_form.lion_id_str.data)
    else:
        out_dct = {}

    return render_template("parser.html", out_dct=out_dct, in_form=in_form)


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
        return jsonify({"code": 1001, "errmsg": "Upload failure!"})


@epilion_blueprint.route("/api/download/<filename>", methods=["GET"])
def download(filename):
    if request.method == "GET":
        if os.path.isfile(os.path.join(app_cfg_dct["ABS_DOWNLOAD_PATH"], filename)):
            return send_from_directory(
                app_cfg_dct["ABS_DOWNLOAD_PATH"], filename, as_attachment=True
            )
        abort(404)


app.register_blueprint(epilion_blueprint)

if __name__ == "__main__":
    app.run()
