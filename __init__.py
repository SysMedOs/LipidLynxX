# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
import json

from flask import Flask
from flask import render_template, redirect, url_for


from epilion.config import DevConfig
from epilion.controllers import Api
from epilion import epilion_blueprint
from epilion.models.DefaultParams import rgx_blank
from epilion.models.Forms import ConverterInputForm
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
    convert_in_form = ConverterInputForm()
    submitted = 0
    if convert_in_form.validate_on_submit():
        usr_abbr_lst = convert_in_form.input_id_str.data.strip("").split("\n")
        usr_abbr_lst = [s for s in usr_abbr_lst if s]

        out_dct = {}
        bad_input_lst = []

        for abbr in usr_abbr_lst:
            converted_info = Api.convert(input_abbreviation=abbr).data
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


app.register_blueprint(epilion_blueprint)

if __name__ == "__main__":
    app.run()
