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

from datetime import datetime
import json
from typing import Union

from fastapi import APIRouter, Request, File, Form, UploadFile
from fastapi.responses import StreamingResponse
from fastapi.templating import Jinja2Templates
import pandas as pd

from lynx.config import api_version, lynx_version
from lynx.models.api_models import (
    FileType,
    StyleType,
    InputListData,
)
import lynx.routers.api as api
from lynx.utils.file_handler import (
    create_converter_output,
    create_equalizer_output, table2html,
)

router = APIRouter()
templates = Jinja2Templates(directory="lynx/templates")


# TemplateResponse from jinja2, ignored in API docs page
@router.get(f"/", include_in_schema=False)
async def home(request: Request):
    return templates.TemplateResponse(
        "home.html",
        {"request": request, "lynx_version": lynx_version, "api_version": api_version},
    )


@router.get("/converter/", include_in_schema=False)
async def converter(request: Request):
    return templates.TemplateResponse(
        "fast_converter.html", {"request": request, "out_dct": {}}
    )


@router.post("/converter/text/", include_in_schema=False)
async def converter_text(
    request: Request,
    lipid_names: str = Form(...),
    export_level: str = Form(...),
    export_style: StyleType = Form(...),
    file_type: FileType = Form(...),
):
    names = lipid_names.splitlines()
    input_data = InputListData(lipid_names=names)
    if export_style == StyleType.lipidlynxx:
        export_style = "LipidLynxX"
    elif export_style == StyleType.comp_db:
        export_style = "COMP_DB"
        export_level = "B1"
    else:
        export_style = "LipidLynxX"
    if file_type == FileType.xlsx:
        file_type = "xlsx"
    elif file_type == FileType.csv:
        file_type = "csv"
    else:
        file_type = "xlsx"
    converted_data = await api.convert_list(
        export_style, input_data, level=export_level
    )
    converted_html, not_converted_html = table2html(converted_data)
    output_name = f"LipidLynxX-Converter_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{file_type}"
    # d = get_download_file(converted_data.data, file_name=output_name, file_type=file_type)
    data_json = json.dumps(converted_data.dict().get("data"))
    print(data_json)
    return templates.TemplateResponse(
        "fast_converter.html",
        {
            "request": request,
            "export_level": export_level,
            "export_style": export_style,
            "output_file_name": output_name,
            "output_file_type": file_type,
            "output_file_data": data_json,
            "converted_html": converted_html,
            "not_converted_html": not_converted_html,
        },
    )


@router.post("/converter/file/", include_in_schema=False)
async def converter_file(
    request: Request,
    file_obj: UploadFile = File(...),
    export_level: str = Form(...),
    export_style: str = Form(...),
    file_type: str = Form(...),
):
    df = pd.read_csv(file_obj.file)
    print(df)
    file_name = "Uploaded"
    return templates.TemplateResponse(
        "fast_converter.html",
        {
            "request": request,
            "file_obj": file_obj,
            "export_level": export_level,
            "export_style": export_style,
            "file_type": file_type,
            "file_name": file_name,
        },
    )


@router.get("/downloads/{data}/{file_type}/{file_name}", name="get_download_file", include_in_schema=False)
def get_download_file(data: str, file_type: str, file_name: str):
    data = json.loads(data)
    if isinstance(data, dict):
        if file_name.startswith("LipidLynxX-Convert"):
            excel_io = create_converter_output(data, file_type=file_type)
        elif file_name.startswith("LipidLynxX-Equal"):
            excel_io = create_equalizer_output(data)
        else:
            try:
                excel_io = create_converter_output(data)
            except Exception as e:
                try:
                    excel_io = create_converter_output(data)
                except Exception as e:
                    raise ValueError()
    else:
        excel_io = None
    if file_name and excel_io:
        return StreamingResponse(
            excel_io, media_type="application/vnd.ms-excel"
        )
