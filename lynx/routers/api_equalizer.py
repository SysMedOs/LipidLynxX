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
from typing import Optional, Union, List

from fastapi import APIRouter, HTTPException, status

from lynx.controllers.equalizer import Equalizer
from lynx.models.api_models import (
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
from lynx.mq.client import equalizer_client
from lynx.utils.job_manager import create_job_token

router = APIRouter()

default_levels = LevelsData(levels=["B1", "D1"])


# POST APIs
@router.post(
    "/single-level/",
    response_model=EqualizerExportData,
    status_code=status.HTTP_201_CREATED,
)
async def equalize_single_level(
    data: InputDictData, level: Optional[str] = "B1"
) -> EqualizerExportData:
    """
    Equalize a dict of lipid names into supported levels and export to supported style
    """
    equalizer = Equalizer(data.data, level=level)
    equalizer_data = equalizer.cross_match()
    return equalizer_data


@router.post(
    "/multiple-levels/",
    response_model=EqualizerExportData,
    status_code=status.HTTP_201_CREATED,
)
async def equalize_multiple_levels(
    data: InputDictData, levels: Optional[LevelsData] = default_levels
) -> EqualizerExportData:
    """
    Equalize a dict of lipid names into supported levels and export to supported style
    """
    # print(levels.levels)
    equalizer = Equalizer(data.data, level=levels.levels)
    equalizer_data = equalizer.cross_match()
    return equalizer_data


@router.post("/dict/", response_model=JobStatus, status_code=status.HTTP_201_CREATED)
async def create_equalize_job(
    data: InputDictData, levels: Union[str, List[str]] = "B1", file_type: str = "xlsx",
):
    """
    """
    token = create_job_token(JobType(job="convert"))
    job_execute_data = {
        "data": data.data,
        "params": {"levels": levels, "file_type": file_type},
    }
    Process(target=equalizer_client, args=(token, job_execute_data, "dict"),).start()

    job_status_data = {
        "data": data.dict(),
        "params": {"levels": levels, "file_type": file_type},
    }
    job_status = "created"
    job_info = JobStatus(token=token, status=job_status, data=job_status_data)

    return job_info
