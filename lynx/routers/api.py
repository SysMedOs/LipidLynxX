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

from fastapi import APIRouter

from lynx.controllers.converter import Converter
from lynx.models.api_models import StyleName, InputListData, ExportListData


router = APIRouter()


# Pure APIs
@router.post("/convert/list/", response_model=ExportListData)
def convert_list(style: StyleName, data: InputListData):
    """
    Convert a list of lipid names into supported abbreviations
    """
    lynx_converter = Converter(style=style)
    print(data)
    x = lynx_converter.convert_list(data.lipid_names)
    return x
