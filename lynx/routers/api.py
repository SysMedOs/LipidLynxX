# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
#
# LipidLynxX is using GPL V3 License
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: a data transfer hub to support integration of large scale lipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import asyncio
import json
import re
from typing import List, Optional, Union
import urllib.parse

from fastapi import APIRouter, HTTPException
from pydantic import parse_obj_as
import requests

from lynx.controllers.linker import get_cross_links
from lynx.controllers.converter import Converter
from lynx.controllers.encoder import Encoder
from lynx.controllers.equalizer import Equalizer
from lynx.controllers.parser import parse_lipid
from lynx.models.api_models import (
    ConverterExportData,
    EqualizerExportData,
    InputDictData,
    InputListData,
    InputStrData,
    LipidNameType,
    LvType,
    LevelsData,
    StyleType,
)
from lynx.models.lipid import LipidType
from lynx.utils.toolbox import get_level

router = APIRouter()

default_levels = LevelsData(levels=["B1", "D1"])


# Get APIs
@router.get("/convert/lipid/")
async def convert_name(
    lipid_name: str = "PLPC",
    style: Optional[str] = "LipidLynxX",
    level: Optional[str] = "MAX",
):
    """
    Convert one lipid name into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = lynx_converter.convert_str(
        input_str=lipid_name, level=get_level(level)
    )
    converted_name = converted_results.output
    if isinstance(converted_name, str) and len(converted_name) > 0:
        pass
    else:
        converted_name = f"Failed to convert: {lipid_name}"

    return converted_name


@router.get("/parse/lipid/")
async def parse_name(lipid_name: str = "PLPC"):
    """
    Parse one lipid name from path parameter
    """
    return parse_lipid(lipid_name)


@router.get("/link/lipid/")
async def link_str(lipid_name: str = "PC(16:0/18:2(9Z,12Z))", export_url: bool = False):
    """
    link one lipid name from data
    """

    linked_ids = await get_cross_links(lipid_name, export_url=export_url)
    if linked_ids:
        return linked_ids
    else:
        raise HTTPException(status_code=500)


# Post APIs
@router.post("/convert/str/")
async def convert_str(
    style: StyleType, data: InputStrData, level: Optional[LvType] = "MAX"
):
    """
    Convert one lipid name into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = lynx_converter.convert_str(
        input_str=data.data, level=get_level(level)
    )
    return converted_results


@router.post("/convert/list/", response_model=ConverterExportData)
async def convert_list(
    style: str, data: InputListData, level: Optional[LvType] = "MAX"
):
    """
    Convert a list of lipid names into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = ConverterExportData(
        data={
            "TextInput": lynx_converter.convert_list(data.data, level=get_level(level))
        }
    )
    return converted_results


@router.post("/convert/dict/", response_model=ConverterExportData)
async def convert_dict(
    style: StyleType, data: InputDictData, level: Optional[LvType] = "MAX"
):
    """
    Convert a dict of lipid names into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = ConverterExportData(
        data=lynx_converter.convert_dict(data.data, level=get_level(level))
    )
    return converted_results


@router.post("/equalize/single-level/", response_model=EqualizerExportData)
async def equalize_single_level(
    data: InputDictData, level: Optional[str] = "B1"
) -> EqualizerExportData:
    """
    Equalize a dict of lipid names into supported levels and export to supported style
    """
    equalizer = Equalizer(data.data, level=level)
    equalizer_data = equalizer.cross_match()
    return equalizer_data


@router.post("/equalize/multiple-levels/", response_model=EqualizerExportData)
async def equalize_multiple_levels(
    data: InputDictData, levels: Optional[LevelsData] = default_levels
) -> EqualizerExportData:
    """
    Equalize a dict of lipid names into supported levels and export to supported style
    """
    print(levels.levels)
    equalizer = Equalizer(data.data, level=levels.levels)
    equalizer_data = equalizer.cross_match()
    return equalizer_data


@router.post("/parse/str/")
async def parse_str(data: InputStrData):
    """
    Parse one lipid name from data
    """

    return parse_lipid(data.data)
