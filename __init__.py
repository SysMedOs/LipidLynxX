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
from flask import abort, jsonify, render_template, redirect, request, url_for
from flask_wtf.csrf import CSRFProtect
from werkzeug.utils import secure_filename

from epilion import app_cfg_dct
from epilion import epilion_blueprint
from epilion.config import DevConfig
from epilion.controllers import Api
from epilion.controllers.FileIO import get_table
from epilion.models.DefaultParams import rgx_blank
from epilion.models.Forms import ConverterTableInputForm
from epilion.models.Forms import ConverterTextInputForm
from epilion.models.Forms import ParserInputForm
from epilion.libLION.LionParser import parse_epilion


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
    convert_in_form = ConverterTextInputForm()
    submitted = 0
    return render_template(
        "converter.html",
        out_dct={},
        bad_in_lst=[],
        in_form=convert_in_form,
        table_form=ConverterTableInputForm(),
        submitted=submitted,
        output_name="",
    )


@epilion_blueprint.route("/converter/text", methods=("GET", "POST"))
def convert_str():
    convert_in_form = ConverterTextInputForm()
    submitted = 0
    if convert_in_form.validate_on_submit():
        usr_abbr_lst = convert_in_form.input_id_str.data.strip("").split("\n")
        usr_abbr_lst = [s for s in usr_abbr_lst if s]

        out_dct = {}
        bad_input_lst = []

        for abbr in usr_abbr_lst:
            converted_info = Api.single_convert(input_abbreviation=abbr).data
            epilion_json = json.loads(converted_info)
            epilion_id = epilion_json["result"]
            if epilion_id:
                out_dct[abbr] = epilion_id
            else:
                if not rgx_blank.match(abbr):
                    bad_input_lst.append(abbr)

        submitted = 1
    else:
        out_dct = {}
        bad_input_lst = []

    if not out_dct:
        bad_input_lst = [""]

    output_name = "test_converter_output.csv"

    return render_template(
        "converter.html",
        out_dct=out_dct,
        bad_in_lst=bad_input_lst,
        in_form=convert_in_form,
        table_form=ConverterTableInputForm(),
        submitted=submitted,
        output_name=output_name,
    )


@epilion_blueprint.route("/converter/table", methods=["POST"])
def convert_table():
    submitted = 0
    out_dct = {}
    table_dct = {}
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
        if masked_file_name:
            masked_path = os.path.join(app_cfg_dct["ABS_UPLOAD_PATH"], masked_file_name)
            usr_file.save(masked_path)
            table_dct = get_table(masked_path)
            submitted = 1
    output_name = "test_converter_output.csv"
    if table_dct:
        out_dct = table_dct
    print({"code": 0, "errmsg": "Upload success.", "result": table_dct})
    return render_template(
        "converter.html",
        out_dct=out_dct,
        bad_in_lst=[],
        in_form=ConverterTextInputForm(),
        table_form=ConverterTableInputForm(),
        submitted=submitted,
        output_name="",
    )


@epilion_blueprint.route("/parser", methods=("GET", "POST"))
def parser():
    in_form = ParserInputForm()
    if in_form.validate_on_submit():
        out_dct = parse_epilion(in_form.lion_id_str.data)
    else:
        out_dct = {}

    return render_template("parser.html", out_dct=out_dct, in_form=in_form)


app.register_blueprint(epilion_blueprint)

if __name__ == "__main__":
    app.run()
