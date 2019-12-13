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
from flask import abort, render_template, redirect, request, url_for
from flask import jsonify, send_from_directory
from flask_restful import Api, Resource, reqparse

from werkzeug.utils import secure_filename

from .config import app_cfg_dct
from .config import blueprint
from .config import DevConfig
from .controllers.encoder import lynx_encode
from .controllers.file_handler import get_table, create_output
from .controllers.parser import parse
from .models.patterns import rgx_blank
from .forms import ConverterTableInputForm
from .forms import ConverterTextInputForm
from .forms import ParserInputForm
from .liblynx.LynxParser import parse_lipidlynx
from .controllers.rest.api import (
    StringConverterAPI,
    DictConverterAPI,
    ListConverterAPI,
    ConverterAPI,
    LevelEqualizerAPI,
    MultiLevelEqualizerAPI,
    EqualizerAPI,
)

app = Flask(__name__)
app.config.from_object(DevConfig)


@app.route("/")
def index():
    return redirect(url_for("lipidlynx.home"))


@blueprint.route("/")
def home():
    return render_template("home.html")


@blueprint.route("/converter", methods=("GET", "POST"))
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


@blueprint.route("/converter/text", methods=("GET", "POST"))
def convert_str():
    convert_in_form = ConverterTextInputForm()
    submitted = 0
    out_df_dct = {"INPUT": [], "OUTPUT": [], "PAIR": [], "NOT_CONVERTED": []}
    if convert_in_form.validate_on_submit():
        usr_abbr_lst = convert_in_form.input_id_str.data.strip("").split("\n")
        usr_abbr_lst = [s for s in usr_abbr_lst if s]

        out_dct = {}
        bad_input_lst = []

        for abbr in usr_abbr_lst:
            converted_info = converter.convert_name(input_abbreviation=abbr).data
            lipidlynx_json = json.loads(converted_info)
            lipidlynx_id = lipidlynx_json["data"]
            if lipidlynx_id:
                out_dct[abbr] = lipidlynx_id
                out_df_dct["INPUT"].append(abbr)
                out_df_dct["OUTPUT"].append(lipidlynx_id)
                out_df_dct["PAIR"].append([abbr, lipidlynx_id])
            else:
                if not rgx_blank.match(abbr):
                    out_df_dct["NOT_CONVERTED"].append(abbr)

        submitted = 1
    else:
        out_dct = {}
        bad_input_lst = []

    output_name = create_output({"ABBREVIATIONS": out_df_dct})

    return render_template(
        "converter.html",
        out_dct=out_dct,
        bad_in_lst=bad_input_lst,
        in_form=convert_in_form,
        table_form=ConverterTableInputForm(),
        submitted=submitted,
        output_name=output_name,
    )


@blueprint.route("/converter/table", methods=["POST"])
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

    if table_dct:
        print({"code": 0, "msg": "Upload success.", "data": table_dct})
        output_info = converter.convert_json(json.dumps(table_dct)).data
        output_json = json.loads(output_info)
        output_dct = output_json["data"]
        output_name = create_output(output_json["data"])
        for k in output_dct:
            pair_lst = output_dct[k].get("PAIR")
            if pair_lst:
                for p in pair_lst:
                    if p[0] not in out_dct and len(p) == 2:
                        out_dct[p[0]] = p[1]
    else:
        output_name = ""

    return render_template(
        "converter.html",
        out_dct=out_dct,
        bad_in_lst=[],
        in_form=ConverterTextInputForm(),
        table_form=ConverterTableInputForm(),
        submitted=submitted,
        output_name=output_name,
    )


@blueprint.route("/parser", methods=("GET", "POST"))
def parser():
    in_form = ParserInputForm()
    if in_form.validate_on_submit():
        out_dct = parse_lipidlynx(in_form.lion_id_str.data)
    else:
        out_dct = {}

    return render_template("parser.html", out_dct=out_dct, in_form=in_form)


# init rest api to blue print
api = Api(blueprint)
api.add_resource(ConverterAPI, "/api/0.1/converter/")
api.add_resource(StringConverterAPI, "/api/0.1/converter/string/")
api.add_resource(ListConverterAPI, "/api/0.1/converter/list/")
api.add_resource(DictConverterAPI, "/api/0.1/converter/dict/")
api.add_resource(EqualizerAPI, "/api/0.1/equalizer/")
api.add_resource(LevelEqualizerAPI, "/api/0.1/equalizer/level/")
api.add_resource(MultiLevelEqualizerAPI, "/api/0.1/equalizer/levels/")

app.register_blueprint(blueprint)


if __name__ == "__main__":

    app.run(debug=True)
