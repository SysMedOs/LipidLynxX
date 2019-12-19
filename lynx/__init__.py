# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os
import requests
import time
from typing import Union, List, Dict

from flask import Flask, send_file
from flask import abort, render_template, redirect, request, url_for
from flask import jsonify, send_from_directory
from flask_restful import Api, Resource, reqparse
import pandas as pd

from werkzeug.utils import secure_filename

from .config import app_cfg_dct
from .config import blueprint
from .config import DevConfig
from .controllers.encoder import lynx_encode
from .controllers.file_handler import get_table, create_output, create_equalizer_output
from .controllers.parser import parse
from .models.defaults import logger, cfg_info_dct
from .models.patterns import rgx_blank
from .forms import (
    ConverterTableInputForm,
    ConverterTextInputForm,
    ConverterForm,
    ParserInputForm,
    EqualizerInputForm,
)
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
# app.config['SERVER_NAME'] = 'lypidlynx.local'
# init rest api to blue print
api = Api(blueprint)
api.add_resource(ConverterAPI, "/api/0.1/converter/")
api.add_resource(StringConverterAPI, "/api/0.1/converter/string/")
api.add_resource(ListConverterAPI, "/api/0.1/converter/list/")
api.add_resource(DictConverterAPI, "/api/0.1/converter/dict/")
api.add_resource(EqualizerAPI, "/api/0.1/equalizer/")
api.add_resource(LevelEqualizerAPI, "/api/0.1/equalizer/level/")
api.add_resource(MultiLevelEqualizerAPI, "/api/0.1/equalizer/levels/")

base_url = cfg_info_dct.get("base_url", "http://127.0.0.1:5000")


def run_converter(input_form):
    submitted = 0
    bad_input_lst = []
    output_name = ""
    excel_data = None
    usr_abbr_lst = input_form.input_id_str.data.strip("").split("\n")
    usr_abbr_lst = [s for s in usr_abbr_lst if s]
    r_url = f"{base_url}{api.url_for(ListConverterAPI)}"
    logger.info(f"Use API ListConverterAPI: {r_url}")
    r = requests.get(r_url, params={"data": json.dumps(usr_abbr_lst)}).json()
    out_dct = r.get("data", {})
    if out_dct:
        submitted = 1
        bad_input_lst = out_dct.get("skipped", [])
        output_name, excel_data = create_output({"ABBREVIATIONS": out_dct})
        logger.info(f"output_name: {output_name}")

    return submitted, out_dct, bad_input_lst, output_name, excel_data


def run_equalizer(data: dict, level: Union[str, List[str]]):
    submitted = 0
    output_name = ""
    r_url = f"{base_url}{api.url_for(EqualizerAPI)}"
    logger.info(f"Use API EqualizerAPI: {r_url}")
    r = requests.get(
        r_url, params={"data": json.dumps(data), "level": json.dumps(level)}
    ).json()
    out_dct = r.get("data", {})
    if out_dct:
        submitted = 1
        output_name = create_equalizer_output(out_dct)
        logger.info(f"output_name: {output_name}")

    return submitted, output_name


@app.route("/")
def index():
    return redirect(url_for("lynx.home"))


@blueprint.route("/")
def home():
    return render_template("home.html")


@blueprint.route("/converter", methods=["GET", "POST"])
def converter():
    return render_template(
        "converter.html",
        out_dct={},
        bad_in_lst=[],
        in_form=ConverterTextInputForm(),
        table_form=ConverterTableInputForm(),
        submitted=0,
        output_name="",
    )


@blueprint.route("/converter/results/string", methods=["GET", "POST"])
def convert_lipid_string():
    usr_input_form = ConverterTextInputForm()

    if usr_input_form.validate_on_submit():
        submitted, out_dct, bad_input_lst, output_name, excel_data = run_converter(
            usr_input_form
        )
    else:
        submitted, out_dct, bad_input_lst, output_name, excel_data = 0, {}, [], "", None
    return render_template(
        "converter.html",
        out_dct=out_dct,
        bad_in_lst=bad_input_lst,
        in_form=usr_input_form,
        table_form=ConverterTableInputForm(),
        submitted=submitted,
        output_name=output_name,
    )


@blueprint.route("/converter/results/table", methods=["GET", "POST"])
def convert_lipid_table():
    logger.info("triggered")
    submitted = 0
    out_dct = {}
    table_dct = {}
    excel_data = {}
    usr_file = request.files["user_file"]
    logger.info(usr_file)
    if usr_file.filename:
        table_dct = get_table(usr_file)
        submitted = 1
    else:
        return render_template(
            "converter.html",
            out_dct={},
            bad_in_lst=["Can not read file or file type not supported."],
            in_form=ConverterTextInputForm(),
            table_form=ConverterTableInputForm(),
            submitted=1,
            output_name="",
        )
    if table_dct:
        logger.info({"code": 0, "msg": "Upload success.", "data": table_dct})
        r_url = f"{base_url}{api.url_for(DictConverterAPI)}"
        logger.info(f"Use API ListConverterAPI: {r_url}")
        r = requests.get(r_url, params={"data": json.dumps(table_dct)}).json()
        excel_data = r["data"]
        output_name = f"LipidLynx_Output_{int(time.time())}.xlsx"
        for k in excel_data:
            pair_lst = excel_data[k].get("converted")
            if pair_lst:
                for p in pair_lst:
                    if p[0] not in out_dct and len(p) == 2:
                        out_dct[p[0]] = p[1]
        excel_data = json.dumps(excel_data)
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
        excel_data=excel_data,
    )


@blueprint.route("/equalizer", methods=("GET", "POST"))
def equalizer():
    submitted = 0
    output_name = ""
    return render_template(
        "equalizer.html",
        in_form=EqualizerInputForm(),
        submitted=submitted,
        output_name=output_name,
    )


@blueprint.route("/equalizer", methods=("GET", "POST"))
def equalize_lipid():
    submitted = 0
    usr_data = {}
    usr_file = request.files["user_file"]
    usr_input_form = EqualizerInputForm()

    if usr_file.filename:
        usr_file_name = secure_filename(usr_file.filename)
        if usr_file_name.endswith(".csv"):
            usr_data = pd.read_csv(usr_file).to_dict(orient="list")
        elif usr_file_name.endswith(".xlsx") or usr_file_name.endswith(".xls"):
            usr_data = pd.read_excel(usr_file).to_dict(orient="list")
        else:
            abort(400, "File type not supported.")

    if usr_data:
        pre_levels = usr_input_form.input_id_str.data.strip("").split("\r\n")

        pre_levels = [s for s in pre_levels if s]
        p_levels = []
        for x in pre_levels:
            p_levels.extend(x.split("\n"))
        levels = []
        for y in p_levels:
            levels.extend(y.split(","))
        logger.info(f"levels {levels}")
        submitted, output_name = run_equalizer(usr_data, level=levels)
    else:
        output_name = ""

    logger.info(output_name)

    return render_template(
        "equalizer.html",
        in_form=usr_input_form,
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


@blueprint.route("/downloads", methods=["GET", "POST"])
def download():
    filename = request.args.get("filename", "")
    data = request.args.get("data", "")
    logger.info(filename)
    if isinstance(data, dict):
        pass
    elif isinstance(data, str):
        data = json.loads(data)

    excel_io = create_output(data)
    if filename and excel_io:
        return send_file(excel_io, attachment_filename=filename, as_attachment=True)


app.register_blueprint(blueprint)


if __name__ == "__main__":

    app.run(debug=True)
