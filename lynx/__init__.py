# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import requests
from datetime import datetime
from typing import Union, List, Dict

from flask import Flask, send_file
from flask import abort, render_template, redirect, request, url_for
import pandas as pd

from werkzeug.utils import secure_filename

# from .liblynx.LynxParser import parse_lipidlynx

from lynx.config import api_url_info, base_url, blueprint, DevConfig
from lynx.forms import (
    ConverterTableInputForm,
    ConverterTextInputForm,
    ConverterForm,
    ParserInputForm,
    EqualizerInputForm,
)
from lynx.models.defaults import logger, cfg_info_dct, api_version
from lynx.utils.file_readers import get_table, create_output, create_equalizer_output
from lynx.utils.toolbox import keep_string_only

lynx_version = 0.4

app = Flask(__name__)
app.config.from_object(DevConfig)
# app.config['SERVER_NAME'] = 'lynx.local'


def run_converter(data: Union[List[str], Dict[str, List[str]]]):
    out_dct = {}
    if isinstance(data, list):
        r_url = api_url_info.get("convert_list")
        logger.info(f"Use API - ListConverterAPI: {r_url}")
    elif isinstance(data, dict):
        r_url = api_url_info.get("convert_dict")
        logger.info(f"Use API - DictConverterAPI: {r_url}")
    else:
        r_url = api_url_info.get("convert")
        logger.info(f"Use API - ConverterAPI: {r_url}")

    r = requests.get(r_url, params={"data": json.dumps(data)}).json()
    excel_data = r["data"]
    output_name = (
        f"LipidLynxX-Converter_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.xlsx"
    )

    for k in excel_data:
        if isinstance(excel_data[k], dict):
            pair_lst = excel_data[k].get("converted")
            if pair_lst:
                for p in pair_lst:
                    if p[0] not in out_dct and len(p) == 2:
                        out_dct[p[0]] = p[1]
        elif isinstance(excel_data[k], list) and k == "converted":
            pair_lst = excel_data.get("converted", [])
            if pair_lst:
                for p in pair_lst:
                    if p[0] not in out_dct and len(p) == 2:
                        out_dct[p[0]] = p[1]
    excel_data = json.dumps(excel_data)

    return out_dct, output_name, excel_data


def run_equalizer(data: dict, level: Union[str, List[str]]):
    submitted = 0
    output_name = (
        f"LipidLynxX-Equalizer_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.xlsx"
    )
    r_url = api_url_info.get("equalize")
    logger.info(f"Use API - EqualizerAPI: {r_url}")
    r = requests.get(
        r_url, params={"data": json.dumps(data), "level": json.dumps(level)}
    ).json()
    excel_data = r.get("data", {})

    if excel_data:
        submitted = 1
        excel_data = json.dumps(r.get("data", {}))
        logger.info("Equalizer")
        logger.info(excel_data)
    else:
        excel_data = None

    return submitted, output_name, excel_data


@app.route("/")
def index():
    return redirect(url_for("lynx.home"))


@blueprint.route("/")
def home():
    return render_template(
        "home.html", lynx_version=lynx_version, api_version=api_version
    )


@blueprint.route("/converter", methods=["GET", "POST"])
def converter():
    return render_template(
        "converter.html",
        out_dct={},
        in_form=ConverterTextInputForm(),
        table_form=ConverterTableInputForm(),
        submitted=0,
        output_name="",
        alerts=[],
    )


@blueprint.route("/converter/results/string", methods=["GET", "POST"])
def convert_lipid_string():
    usr_input_form = ConverterTextInputForm()

    if usr_input_form.validate_on_submit():
        usr_data = usr_input_form.input_id_str.data.strip("").split("\n")
        usr_data = [s for s in usr_data if s]
        out_dct, output_name, excel_data = run_converter(usr_data)
        submitted = 1
    else:
        submitted, out_dct, bad_input_lst, output_name, excel_data = 0, {}, [], "", None
    return render_template(
        "converter.html",
        out_dct=out_dct,
        in_form=usr_input_form,
        table_form=ConverterTableInputForm(),
        submitted=submitted,
        output_name=output_name,
        excel_data=excel_data,
        alerts=[],
    )


@blueprint.route("/converter/results/table", methods=["GET", "POST"])
def convert_lipid_table():
    logger.info("triggered")
    out_dct = {}
    excel_data = {}
    usr_file = request.files["user_file"]
    logger.info(usr_file)
    usr_table_form = ConverterTableInputForm()
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
            excel_data=excel_data,
            alerts=[],
        )
    table_dct = keep_string_only(table_dct)
    if table_dct:
        logger.info({"code": 0, "msg": "Upload success.", "data": table_dct})
        out_dct, output_name, excel_data = run_converter(table_dct)
        submitted = 1
    else:
        output_name = ""

    return render_template(
        "converter.html",
        out_dct=out_dct,
        in_form=ConverterTextInputForm(),
        table_form=usr_table_form,
        submitted=submitted,
        output_name=output_name,
        excel_data=excel_data,
        alerts=[],
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
        excel_data=None,
        alert_lst=[],
    )


@blueprint.route("/equalizer/results", methods=("GET", "POST"))
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
    else:
        return render_template(
            "equalizer.html",
            in_form=EqualizerInputForm(),
            submitted=1,
            output_name="",
            alerts=["Can not read file or file type not supported."],
        )

    if usr_data:
        pre_levels = usr_input_form.input_id_str.data.strip("").split("\r\n")

        pre_levels = [s for s in pre_levels if s]
        p_levels = []
        for x in pre_levels:
            p_levels.extend(x.split("\n"))
        levels = []
        for y in p_levels:
            levels.extend(y.split(","))
        levels = [lv.strip(" ") for lv in levels]
        logger.info(f"levels {levels}")
        submitted, output_name, excel_data = run_equalizer(usr_data, level=levels)
    else:
        output_name = ""
        excel_data = None

    return render_template(
        "equalizer.html",
        in_form=usr_input_form,
        submitted=submitted,
        output_name=output_name,
        alerts=[],
        excel_data=excel_data,
    )


# @blueprint.route("/parser", methods=("GET", "POST"))
# def parser():
#     in_form = ParserInputForm()
#     if in_form.validate_on_submit():
#         out_dct = parse_lipidlynx(in_form.lion_id_str.data)
#     else:
#         out_dct = {}
#
#     return render_template("parser.html", out_dct=out_dct, in_form=in_form)


@blueprint.route("/downloads", methods=["GET", "POST"])
def download():
    filename = request.args.get("filename", "")
    data = request.args.get("data", "")
    logger.info(filename)
    if isinstance(data, dict):
        pass
    elif isinstance(data, str):
        data = json.loads(data)
    if isinstance(data, dict):
        if filename.startswith("LipidLynxX-Convert"):
            excel_io = create_output(data)
        elif filename.startswith("LipidLynxX-Equal"):
            excel_io = create_equalizer_output(data)
        else:
            try:
                excel_io = create_output(data)
            except Exception as e:
                try:
                    excel_io = create_output(data)
                except Exception as e:
                    raise ValueError()
    else:
        excel_io = None
    if filename and excel_io:
        return send_file(
            excel_io,
            mimetype="application/vnd.ms-excel",
            attachment_filename=filename,
            as_attachment=True,
        )


@blueprint.route("/user_guide")
def user_guide():
    return render_template(
        "user_guide.html", lynx_version=lynx_version, api_version=api_version
    )


@blueprint.route("/nomenclature")
def nomenclature():
    return render_template(
        "nomenclature.html"
    )


@blueprint.route("/levels")
def levels():
    return render_template(
        "levels.html"
    )


@blueprint.route("/about")
def about():
    return render_template(
        "about.html", lynx_version=lynx_version, api_version=api_version
    )


app.register_blueprint(blueprint)


if __name__ == "__main__":

    app.run(debug=True)
