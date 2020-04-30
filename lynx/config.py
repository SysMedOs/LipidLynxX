# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# LipidLynxX is Dual-licensed
#   For academic and non-commercial use: GPLv2 License:
#   For commercial use: please contact the SysMedOs team by email.
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: lipid annotations converter for large scale lipidomics and epilipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os

from flask import Blueprint
from flask_restful import Api

from lynx.controllers.rest.api import (
    StringConverterAPI,
    DictConverterAPI,
    ListConverterAPI,
    ConverterAPI,
    LevelEqualizerAPI,
    MultiLevelEqualizerAPI,
    EqualizerAPI,
)
from lynx.models.defaults import logger, cfg_info_dct, api_version


class Config(object):
    SECRET_KEY = os.urandom(24).hex()


class ProdConfig(Config):
    pass


class DevConfig(Config):
    DEBUG = True


base_url = cfg_info_dct.get("base_url", "http://127.0.0.1:5000")
app_name = "lynx"
app_url_prefix = "/lynx"

blueprint = Blueprint(
    app_name,
    __name__,
    template_folder=r"templates",
    static_folder=r"static",
    url_prefix=app_url_prefix,
)

# init rest api to blue print
api = Api(blueprint)
api_cfg_info = {
    "convert": (ConverterAPI, f"/api/{api_version}/converter/"),
    "convert_str": (StringConverterAPI, f"/api/{api_version}/converter/str/"),
    "convert_list": (ListConverterAPI, f"/api/{api_version}/converter/list/"),
    "convert_dict": (DictConverterAPI, f"/api/{api_version}/converter/dict/"),
    "equalize": (EqualizerAPI, f"/api/{api_version}/equalizer/"),
    "equalize_to_level": (LevelEqualizerAPI, f"/api/{api_version}/equalizer/level/"),
    "equalize_to_levels": (
        MultiLevelEqualizerAPI,
        f"/api/{api_version}/equalizer/levels/",
    ),
}
api_url_info = {}
for api_name in api_cfg_info:
    api.add_resource(api_cfg_info[api_name][0], api_cfg_info[api_name][1])
    api_url_info[api_name] = f"{base_url}{app_url_prefix}{api_cfg_info[api_name][1]}"

logger.info("Lynx API config loaded...")

# app_cfg_dct = {
#     "ABS_BASE_PATH": os.path.abspath(os.path.dirname(__file__)),
#     "UPLOAD_FOLDER": "uploads",
#     "DOWNLOAD_FOLDER": "downloads",
# }
#
# app_cfg_dct["ABS_UPLOAD_PATH"] = os.path.join(
#     app_cfg_dct["ABS_BASE_PATH"], app_cfg_dct["UPLOAD_FOLDER"]
# )
# app_cfg_dct["ABS_DOWNLOAD_PATH"] = os.path.join(
#     app_cfg_dct["ABS_BASE_PATH"], app_cfg_dct["DOWNLOAD_FOLDER"]
# )
