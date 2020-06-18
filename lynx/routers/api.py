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
import json
from typing import List, Optional, Union

from fastapi import APIRouter

from lynx.controllers.converter import Converter
from lynx.controllers.encoder import Encoder
from lynx.controllers.equalizer import Equalizer
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
from lynx.utils.toolbox import get_level

router = APIRouter()

default_levels = LevelsData(levels=["B1", "D1"])


# Get APIs
@router.get("/convert/lipid/")
async def convert_name(
    lipid_name: str, style: Optional[str] = "LipidLynxX", level: Optional[str] = "MAX"
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
    # extractor = Decoder()
    # parsed = extractor.extract(data)
    lynx_gen = Encoder()
    parsed = lynx_gen.export_all_levels(lipid_name)

    return parsed


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
            "TextInput": lynx_converter.convert_list(
                data.data, level=get_level(level)
            )
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
async def parse_str(data: LipidNameType = "PLPC"):
    """
    Parse one lipid name from data
    """
    # extractor = Decoder()
    # parsed = extractor.extract(data)
    lynx_gen = Encoder()
    parsed = lynx_gen.export_all_levels(data)

    return parsed
