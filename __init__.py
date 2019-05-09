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

from epilion.config import DevConfig

epilion_blueprint = Blueprint(
    'epilion',
    __name__,
    template_folder=r'epilion/templates',
    static_folder=r'epilion/static',
    url_prefix='/epilion'
)

app = Flask(__name__)
app.config.from_object(DevConfig)


print(os.getcwd())


@app.route('/')
def index():
    return redirect(url_for('epilion.home'))


@epilion_blueprint.route('/')
@epilion_blueprint.route('/<int:page>')
def home():
    return render_template('home.html')


app.register_blueprint(epilion_blueprint)

if __name__ == '__main__':
    app.run()
