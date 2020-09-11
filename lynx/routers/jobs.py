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
import time
from typing import Union, Optional

from fastapi import APIRouter, HTTPException, status
from starlette import status
from starlette.background import BackgroundTasks
from starlette.concurrency import run_in_threadpool

from lynx import Converter
from lynx.models.api_models import (
    JobStatus,
    InputDictData,
    InputListData,
    InputStrData,
    StyleType,
    LvType,
    ConverterExportData,
    JobType,
)
from lynx.routers.api import link_one_lipid
from lynx.utils.job_manager import (
    is_job_finished,
    get_job_output,
    save_session,
    create_job_token,
)
from lynx.utils.toolbox import get_level
from lynx.tasks.client import converter_client


router = APIRouter()


@router.get(
    "/converter/{token}",
    response_model=JobStatus,
    status_code=status.HTTP_202_ACCEPTED,
)
async def check_converter_job(token: str,):
    """"""
    if is_job_finished(token):
        job_data = get_job_output(token)
        job_status = "finished"
    else:
        job_data = {}
        job_status = "working"
    job_info = JobStatus(token=token, status=job_status, data=job_data)
    return job_info


@router.get(
    "/linker/{token}", response_model=JobStatus, status_code=status.HTTP_202_ACCEPTED,
)
async def check_linker_job(token: str,):
    """"""
    if is_job_finished(token):
        job_data = get_job_output(token)
        job_status = "finished"
    else:
        job_data = {}
        job_status = "working"
    job_info = JobStatus(token=token, status=job_status, data=job_data)
    return job_info


async def run_converter(
    token: str,
    data: Union[InputDictData, InputListData, InputStrData],
    style: StyleType,
    level: Optional[LvType] = "MAX",
):
    """
    link a list of lipids to related resources from posted lipid name list
    """

    lynx_converter = Converter(style=style)
    if isinstance(data, InputStrData):
        converted_results = lynx_converter.convert_str(
            input_str=data.data, level=get_level(level)
        )  # ConverterExportData
        save_session(token, data=converted_results.dict())
    elif isinstance(data, InputListData):
        results = await run_in_threadpool(
            lynx_converter.convert_list, data.data, level=get_level(level)
        )
        converted_results = ConverterExportData(data={"TextInput": results})
        # converted_results = ConverterExportData(
        #     data={
        #         "TextInput": lynx_converter.convert_list(
        #             data.data, level=get_level(level)
        #         )
        #     }
        # )
        save_session(token, data=converted_results.dict())
    elif isinstance(data, InputDictData):
        converted_results = ConverterExportData(
            data=lynx_converter.convert_dict(data.data, level=get_level(level))
        )
        save_session(token, data=converted_results.dict())


@router.post("/convert/", response_model=JobStatus, status_code=status.HTTP_201_CREATED)
async def create_convert_job(
    background_tasks: BackgroundTasks,
    data: Union[InputDictData, InputListData, InputStrData],
    style: StyleType,
    level: Optional[LvType] = "MAX",
):
    """"""
    token = create_job_token(JobType(job="convert"))

    client_data = {
        "names": data.data,
        "export_style": style,
        "export_level": level,
    }
    Process(
        target=converter_client,
        args=(
            token,
            client_data
        ),
    ).start()

    job_status_data = {
        "data": data.dict(),
        "style": style,
        "level": level,
    }
    job_execute_data = {
        "data": data,
        "style": style,
        "level": level,
    }
    job_status = "created"
    job_info = JobStatus(token=token, status=job_status, data=job_status_data)
    background_tasks.add_task(run_converter, token, **job_execute_data)
    return job_info


async def run_linker(
    token: str, lipid_names: list, export_url: bool = False, export_names: bool = True,
):
    """
    link a list of lipids to related resources from posted lipid name list
    """
    linked_info = {}
    for lipid_name in lipid_names:
        linked_info[lipid_name] = await link_one_lipid(
            lipid_name, export_url, export_names
        )
    save_session(token, data=linked_info)


@router.post("/link/", response_model=JobStatus, status_code=status.HTTP_201_CREATED)
async def create_linker_job(
    lipid_names: list,
    background_tasks: BackgroundTasks,
    export_url: bool = False,
    export_names: bool = True,
):
    """"""
    token = create_job_token(JobType(job="link"))
    job_data = {
        "data": lipid_names,
        "export_url": export_url,
        "export_names": export_names,
    }
    job_status = "created"
    job_info = JobStatus(token=token, status=job_status, data=job_data)
    background_tasks.add_task(
        run_linker, token, lipid_names, export_url=export_url, export_names=export_names
    )
    return job_info
