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

from fastapi import APIRouter, HTTPException, status

from lynx.controllers.linker import get_cross_links, get_lmsd_name, get_swiss_name
from lynx.models.api_models import JobStatus, JobType, LevelsData, InputListData
from lynx.mq.client import linker_client
from lynx.routers.api_converter import convert_lipid
from lynx.utils.job_manager import (
    is_job_finished,
    get_job_output,
    save_job,
    create_job_token,
)

router = APIRouter()

default_levels = LevelsData(levels=["B1", "D1"])


# Reusable functions
async def link_one_lipid(
    lipid_name: str = "PC(16:0/18:2(9Z,12Z))",
    export_url: bool = False,
    export_names: bool = True,
):
    if lipid_name:
        if re.match(r"^LM\w\w\d{8}$", lipid_name, re.IGNORECASE):
            safe_lipid_name = await get_lmsd_name(lipid_name)
        elif re.match(r"^SLM:\d{9}$", lipid_name, re.IGNORECASE):
            safe_lipid_name = await get_swiss_name(lipid_name)
        else:
            safe_lipid_name = lipid_name
    else:
        raise HTTPException(status_code=404)
    search_name = await convert_lipid(
        safe_lipid_name, level="MAX", style="BracketsShorthand"
    )
    if (
        search_name.startswith("SM(")
        and not search_name.startswith("SM(d")
        and not search_name.startswith("SM(m")
        and not search_name.startswith("SM(t")
    ):
        if ";O" in search_name:
            search_name = re.sub(r"\(", "(m", search_name)
            search_name = re.sub(r";O", "", search_name)
        elif ";1" in search_name:
            search_name = re.sub(r"\(", "(m", search_name)
            search_name = re.sub(r";1", "", search_name)
        elif ";2" in search_name:
            search_name = re.sub(r"\(", "(d", search_name)
            search_name = re.sub(r";2", "", search_name)
        elif ";3" in search_name:
            search_name = re.sub(r"\(", "(t", search_name)
            search_name = re.sub(r";3", "", search_name)
    elif (
        search_name.startswith("Cer(")
        and not search_name.startswith("Cer(d")
        and not search_name.startswith("Cer(m")
        and not search_name.startswith("Cer(t")
    ):
        if ";O" in search_name:
            search_name = re.sub(r"\(", "(m", search_name)
            search_name = re.sub(r";O", "", search_name)
        elif ";1" in search_name:
            search_name = re.sub(r"\(", "(m", search_name)
            search_name = re.sub(r";1", "", search_name)
        elif ";2" in search_name:
            search_name = re.sub(r"\(", "(d", search_name)
            search_name = re.sub(r";2", "", search_name)
        elif ";3" in search_name:
            search_name = re.sub(r"\(", "(t", search_name)
            search_name = re.sub(r";3", "", search_name)
    elif (
        search_name.startswith("HexCer(")
        and not search_name.startswith("HexCer(d")
        and not search_name.startswith("HexCer(m")
        and not search_name.startswith("HexCer(t")
    ):
        if ";O" in search_name:
            search_name = re.sub(r"\(", "(m", search_name)
            search_name = re.sub(r";O", "", search_name)
        elif ";1" in search_name:
            search_name = re.sub(r"\(", "(m", search_name)
            search_name = re.sub(r";1", "", search_name)
        elif ";2" in search_name:
            search_name = re.sub(r"\(", "(d", search_name)
            search_name = re.sub(r";2", "", search_name)
        elif ";3" in search_name:
            search_name = re.sub(r"\(", "(t", search_name)
            search_name = re.sub(r";3", "", search_name)
    elif (
        search_name.startswith("Hex2Cer(")
        and not search_name.startswith("Hex2Cer(d")
        and not search_name.startswith("Hex2Cer(m")
        and not search_name.startswith("Hex2Cer(t")
    ):
        if ";O" in search_name:
            search_name = re.sub(r"\(", "(m", search_name)
            search_name = re.sub(r";O", "", search_name)
        elif ";1" in search_name:
            search_name = re.sub(r"\(", "(m", search_name)
            search_name = re.sub(r";1", "", search_name)
        elif ";2" in search_name:
            search_name = re.sub(r"\(", "(d", search_name)
            search_name = re.sub(r";2", "", search_name)
        elif ";3" in search_name:
            search_name = re.sub(r"\(", "(t", search_name)
            search_name = re.sub(r";3", "", search_name)
    print("search_name", search_name)
    resource_data = await get_cross_links(search_name, export_url=export_url)
    if resource_data:
        if export_names:
            render_data_dct = {
                "lipid_name": lipid_name,
                "search_name": search_name,
                "shorthand_name": await convert_lipid(
                    safe_lipid_name, level="MAX", style="ShorthandNotation"
                ),
                "lynx_name": await convert_lipid(
                    safe_lipid_name, level="MAX", style="LipidLynxX"
                ),
                "biopan_name": await convert_lipid(
                    safe_lipid_name, level="B2", style="BioPAN"
                ),
                "resource_data": resource_data,
            }
            return render_data_dct
        else:
            return resource_data
    else:
        raise HTTPException(status_code=500)


@router.get("/lipid/")
async def get_link_lipid(
    lipid_name: str = "PC(16:0/18:2(9Z,12Z))",
    export_url: bool = False,
    export_names: bool = True,
) -> dict:
    """
    link one lipid  to related resources from given lipid name
    """
    return await link_one_lipid(lipid_name, export_url, export_names)


@router.get(
    "/jobs/{token}", response_model=JobStatus, status_code=status.HTTP_202_ACCEPTED,
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


# Post APIs
@router.post("/str/", status_code=status.HTTP_201_CREATED)
async def link_str(
    lipid_name: str = "PC(16:0/18:2(9Z,12Z))",
    export_url: bool = False,
    export_names: bool = True,
) -> dict:
    """
    link one lipid to related resources from posted lipid name
    """
    return await link_one_lipid(lipid_name, export_url, export_names)


@router.post("/list/", status_code=status.HTTP_201_CREATED)
async def link_list(
    lipid_names: list, export_url: bool = False, export_names: bool = True,
) -> dict:
    """
    link a list of lipids to related resources from posted lipid name list
    """
    linked_info = {}
    for lipid_name in lipid_names:
        linked_info[lipid_name] = await link_one_lipid(
            lipid_name, export_url, export_names
        )
    return linked_info


@router.post("/dict/", status_code=status.HTTP_201_CREATED)
async def link_dict(
    lipid_names: dict, export_url: bool = False, export_names: bool = True
) -> dict:
    """
    link a list of lipids to related resources from posted lipid name list
    """
    linked_info = {}
    for col in lipid_names:
        lipid_col = lipid_names[col]
        linked_col_info = {}
        for lipid_name in lipid_col:
            linked_col_info[lipid_name] = await link_one_lipid(
                lipid_name, export_url, export_names
            )
        linked_info[col] = linked_col_info
    return linked_info


@router.post("/list/", response_model=JobStatus, status_code=status.HTTP_201_CREATED)
async def create_link_list_job(
    data: InputListData,
    export_url: bool = True,
    export_names: bool = True,
    file_type: str = "xlsx",
):
    """"""
    token = create_job_token(JobType(job="link"))
    job_execute_data = {
        "data": data.data,
        "params": {
            "export_url": export_url,
            "export_names": export_names,
            "file_type": file_type,
        },
    }
    Process(target=linker_client, args=(token, job_execute_data, "list"),).start()
    job_status = "created"
    job_info = JobStatus(token=token, status=job_status, data=job_execute_data)

    return job_info
