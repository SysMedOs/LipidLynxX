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
    resource_data = await get_cross_links(search_name, export_url=export_url)
    if resource_data:
        if export_names:
            render_data_dct = {
                "lipid_name": lipid_name,
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
async def link_lipid(
    lipid_name: str = "PC(16:0/18:2(9Z,12Z))",
    export_url: bool = False,
    export_names: bool = True,
) -> dict:
    """
    link one lipid  to related resources from given lipid name
    """
    return await link_one_lipid(lipid_name, export_url, export_names)


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
    lipid_names: list,
    export_url: bool = False,
    export_names: bool = True,
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
