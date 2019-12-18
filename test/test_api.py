# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import requests

from lynx import api, app
from lynx.controllers.rest.api import (
    StringConverterAPI,
    DictConverterAPI,
    ListConverterAPI,
    ConverterAPI,
    LevelEqualizerAPI,
    MultiLevelEqualizerAPI,
    EqualizerAPI,
)

# app.run(debug=True)

# print(api.url_for(StringConverterAPI))
r = requests.get("localhost:5000/api/0.1/converter/", data="PC 18:0_18:2")
print(r)
