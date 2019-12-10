# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os

from flask import Blueprint

from .controllers.encoder import lynx_encode
from .controllers.parser import parse


class Config(object):
    SECRET_KEY = os.urandom(24).hex()


class ProdConfig(Config):
    pass


class DevConfig(Config):
    DEBUG = True


app_cfg_dct = {
    "ABS_BASE_PATH": os.path.abspath(os.path.dirname(__file__)),
    "UPLOAD_FOLDER": "uploads",
    "DOWNLOAD_FOLDER": "downloads",
}

app_cfg_dct["ABS_UPLOAD_PATH"] = os.path.join(
    app_cfg_dct["ABS_BASE_PATH"], app_cfg_dct["UPLOAD_FOLDER"]
)
app_cfg_dct["ABS_DOWNLOAD_PATH"] = os.path.join(
    app_cfg_dct["ABS_BASE_PATH"], app_cfg_dct["DOWNLOAD_FOLDER"]
)

lipidlynx_blueprint = Blueprint(
    "lipidlynx",
    __name__,
    template_folder=r"templates",
    static_folder=r"static",
    url_prefix="/lipidlynx",
)
