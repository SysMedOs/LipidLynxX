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

from multiprocessing import Process
import re
from typing import Optional

from fastapi import APIRouter, HTTPException, status

from lynx.controllers.linker import get_cross_links, get_lmsd_name, get_swiss_name
from lynx.controllers.converter import Converter
from lynx.controllers.equalizer import Equalizer
from lynx.controllers.parser import parse_lipid
from lynx.models.api_models import (
    ConverterExportData,
    EqualizerExportData,
    InputDictData,
    InputListData,
    InputStrData,
    JobStatus,
    JobType,
    LvType,
    LevelsData,
    StyleType,
)
from lynx.models.defaults import (
    default_temp_folder,
    default_temp_max_days,
    default_temp_max_files,
)
from lynx.utils.log import app_logger
from lynx.utils.toolbox import get_level
from lynx.utils.file_handler import clean_temp_folder
from lynx.tasks.client import converter_client
from lynx.utils.job_manager import create_job_token


router = APIRouter()

default_levels = LevelsData(levels=["B1", "D1"])

# Get APIs
@router.get("/lipid/")
async def convert_lipid(
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


# Post APIs
@router.post("/str/", response_model=JobStatus, status_code=status.HTTP_201_CREATED)
async def convert_str(
    style: StyleType,
    data: InputStrData,
    level: Optional[LvType] = "MAX",
):
    """
    Convert one lipid name into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = lynx_converter.convert_str(
        input_str=data.data, level=get_level(level)
    )
    return converted_results


@router.post(
    "/list/",
    response_model=ConverterExportData,
    status_code=status.HTTP_201_CREATED,
)
def convert_list(data: InputListData, style: str, level: Optional[LvType] = "MAX"):
    """
    Convert a list of lipid names into supported levels and export to supported style
    """
    token = create_job_token(JobType(job="convert"))
    client_data = {}
    Process(
        target=converter_client,
        args=(token, client_data),
    ).start()
    # lynx_converter = Converter(style=style)
    # converted_results = ConverterExportData(
    #     data={
    #         "TextInput": lynx_converter.convert_list(data.data, level=get_level(level))
    #     }
    # )
    # return converted_results


@router.post(
    "/dict/",
    response_model=ConverterExportData,
    status_code=status.HTTP_201_CREATED,
)
async def convert_dict(
    data: InputDictData, style: StyleType, level: Optional[LvType] = "MAX"
):
    """
    Convert a dict of lipid names into supported levels and export to supported style
    """
    lynx_converter = Converter(style=style)
    converted_results = ConverterExportData(
        data=lynx_converter.convert_dict(data.data, level=get_level(level))
    )
    return converted_results
