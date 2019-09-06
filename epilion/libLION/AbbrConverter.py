# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from typing import Tuple

from epilion.controllers.DefaultParams import abbr_cfg_path
from epilion.libLION.Converter import Converter


def web_converter(usr_input: str) -> Tuple[dict, list]:
    converter = Converter(abbr_cfg_path)
    epilion_dct, bad_input_lst = converter.convert_text(usr_input)

    return epilion_dct, bad_input_lst
