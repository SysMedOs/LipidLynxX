# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
import os

from flask import Flask
from flask import Blueprint
from flask import render_template, redirect, url_for
from flask_wtf import Form
from wtforms import StringField, TextAreaField
from wtforms.validators import DataRequired, Length

from epilion.config import DevConfig
from epilion.controllers.AbbrConverter import web_converter
from epilion.controllers.LionParser import parse_epilion


class ConverterInputForm(Form):

    input_id_str = TextAreaField(u'Paste lipid abbreviations here:', validators=[DataRequired()])


class ParserInputForm(Form):

    lion_id_str = StringField(u'Paste epiLION id here:', validators=[DataRequired(), Length(max=255)])


epilion_blueprint = Blueprint(
    'epilion',
    __name__,
    template_folder=r'epilion/templates',
    static_folder=r'epilion/static',
    url_prefix='/epilion'
)

app = Flask(__name__)
app.config.from_object(DevConfig)


@app.route('/')
def index():
    return redirect(url_for('epilion.home'))


@epilion_blueprint.route('/')
def home():
    return render_template('home.html')


@epilion_blueprint.route('/converter', methods=('GET', 'POST'))
def converter():
    conv_in_form = ConverterInputForm()
    if conv_in_form.validate_on_submit():
        out_lst = web_converter(conv_in_form.input_id_str.data)
    else:
        out_lst = ['']

    return render_template('converter.html', out_lst=out_lst, in_form=conv_in_form)


@epilion_blueprint.route('/parser', methods=('GET', 'POST'))
def parser():
    in_form = ParserInputForm()
    if in_form.validate_on_submit():
        out_dct = parse_epilion(in_form.lion_id_str.data)
    else:
        out_dct = {}

    return render_template('parser.html', out_dct=out_dct, in_form=in_form)


app.register_blueprint(epilion_blueprint)

if __name__ == '__main__':
    app.run()
