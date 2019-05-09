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
from wtforms import  StringField, TextAreaField
from wtforms.validators import DataRequired, Length

from epilion.config import DevConfig
from epilion.controllers.AbbrConverter import web_converter


class ConverterInputForm(Form):

    input_id_str = TextAreaField(u'Paste lipid abbreviations here:', validators=[DataRequired()])


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

    # out_str = '\n'.join(out_lst)

    return render_template('converter.html', out_lst=out_lst, in_form=conv_in_form)


app.register_blueprint(epilion_blueprint)

if __name__ == '__main__':
    app.run()
