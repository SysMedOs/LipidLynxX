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
from lynx.utils.temp_file_cleaner import clean_temp_folder
from lynx.mq.client import converter_client
from lynx.utils.job_manager import create_job_token


router = APIRouter()

default_levels = LevelsData(levels=["B1", "D1"])


# Get APIs
@router.get("/lipid/")
async def parse_lipid(lipid_name: str = "PLPC"):
    """
    Parse one lipid name from path parameter
    """
    return parse_lipid(lipid_name)


# Post APIs
@router.post("/str/", status_code=status.HTTP_201_CREATED)
async def parse_str(data: InputStrData):
    """
    Parse one lipid name from data
    """

    return parse_lipid(data.data)
