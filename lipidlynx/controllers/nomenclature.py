# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from collections import Counter
import re
from typing import Dict, List, Tuple, Union

from natsort import natsorted

from lipidlynx.models.defaults import (
    lipid_class_alias_info,
    cv_order_list,
    cv_alias_info,
)
from lipidlynx.controllers.Logger import logger
from lipidlynx.controllers.GeneralFunctions import seg_to_str
from lipidlynx.controllers.Parser import parse, parse_mod


class LynxObject(object):

    def __init__(self, lynx_id: str):

        self.lynx_id = lynx_id

    def structure(self):
        pass
