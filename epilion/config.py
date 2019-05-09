# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de


class Config(object):
    SECRET_KEY = '1472130a246ecd04a61ef528fd5151de'


class ProdConfig(Config):
    pass


class DevConfig(Config):
    DEBUG = True
