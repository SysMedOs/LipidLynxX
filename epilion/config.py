# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from .configurations.key_cfg import secrete_key


class Config(object):
    SECRET_KEY = secrete_key


class ProdConfig(Config):
    pass


class DevConfig(Config):
    DEBUG = True
