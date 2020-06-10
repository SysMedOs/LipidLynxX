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

from typing import Dict, List, Optional, Union

from fastapi import APIRouter

from lynx.controllers.converter import Converter
from lynx.controllers.decoder import Decoder
from lynx.controllers.encoder import Encoder
from lynx.models.api_models import (
    LvType,
    LipidNameType,
    StyleType,
    InputDictData,
    InputListData,
    InputStrData,
    ConverterExportDictData,
    ConverterExportListData,
)


router = APIRouter()


def get_level(lv) -> str:
    if isinstance(lv, LvType) or isinstance(lv, str):
        use_level = lv
    else:
        use_level = "MAX"
    return use_level


# Pure APIs
@router.get("/parse/str/")
def parse_str(data: LipidNameType = "PLPC"):
    """
    Parse one lipid name
    """
    # extractor = Decoder()
    # parsed = extractor.extract(data)
    lynx_gen = Encoder()
    parsed = lynx_gen.export_all_levels(data)

    return parsed


@router.post("/convert/str/")
def convert_str(style: StyleType, data: InputStrData, level: Optional[LvType] = "MAX"):
    """
    Convert one lipid name into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = lynx_converter.convert_str(
        input_str=data.lipid_name, level=get_level(level)
    )
    return converted_results


@router.post("/convert/list/", response_model=ConverterExportListData)
def convert_list(
    style: StyleType, data: InputListData, level: Optional[LvType] = "MAX"
):
    """
    Convert a list of lipid names into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = lynx_converter.convert_list(
        data.lipid_names, level=get_level(level)
    )
    return converted_results


@router.post("/convert/dict/", response_model=ConverterExportDictData)
def convert_dict(
    style: StyleType, data: InputDictData, level: Optional[LvType] = "MAX"
):
    """
    Convert a dict of lipid names into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = {
        "data": lynx_converter.convert_dict(data.data, level=get_level(level))
    }
    return converted_results
