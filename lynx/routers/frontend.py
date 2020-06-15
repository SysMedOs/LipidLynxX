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

import base64
import json
from typing import IO
from tempfile import NamedTemporaryFile

from fastapi import (
    APIRouter,
    File,
    Form,
    Request,
    UploadFile,
)
from fastapi.responses import StreamingResponse
from fastapi.templating import Jinja2Templates


from lynx.config import api_version, lynx_version
from lynx.models.api_models import FileType, StyleType, InputListData, InputDictData
import lynx.routers.api as api
from lynx.utils.file_handler import (
    create_converter_output,
    create_equalizer_output,
    get_table,
    get_file_type,
    get_output_name,
    table2html,
)
from lynx.utils.toolbox import get_levels, get_style_level, get_url_safe_str


# upload file size check
# async def valid_content_length(content_length: int = Header(..., lt=10240_000)):
#     return content_length


router = APIRouter()
templates = Jinja2Templates(directory="lynx/templates")


# TemplateResponse from jinja2, ignored in API docs page
@router.get("/about", include_in_schema=False)
def about(request: Request):
    return templates.TemplateResponse(
        "about.html",
        {"request": request, "lynx_version": lynx_version, "api_version": api_version},
    )


@router.get("/converter/", include_in_schema=False)
async def converter(request: Request):
    return templates.TemplateResponse(
        "converter.html", {"request": request, "out_dct": {}}
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
    export_style, export_level = get_style_level(export_style, export_level)
    converted_data = await api.convert_list(export_style, input_data, export_level)
    converted_html, not_converted_html = table2html(converted_data)
    data_encoded = get_url_safe_str(converted_data.dict().get("data"))
    file_type = get_file_type(file_type)
    output_name = get_output_name("Converter", file_type)
    render_data_dct = {
        "request": request,
        "export_level": export_level,
        "export_style": export_style,
        "output_file_name": output_name,
        "output_file_type": file_type,
        "output_file_data": data_encoded,
        "converted_html": converted_html,
        "not_converted_html": not_converted_html,
    }
    return templates.TemplateResponse("converter.html", render_data_dct)


@router.post("/converter/file/", include_in_schema=False)
async def converter_file(
    request: Request,
    file_obj: UploadFile = File(...),
    export_level: str = Form(...),
    export_style: StyleType = Form(...),
    file_type: FileType = Form(...),
):
    table_info, err_lst = get_table(file_obj, err_lst=[])
    if table_info:
        export_style, export_level = get_style_level(export_style, export_level)
        input_data = InputDictData(data=table_info)
        converted_data = await api.convert_dict(export_style, input_data, export_level)
        converted_html, not_converted_html = table2html(converted_data)
        data_encoded = get_url_safe_str(converted_data.dict().get("data"))
        file_type = get_file_type(file_type)
        output_name = get_output_name("Converter", file_type)
        render_data_dct = {
            "request": request,
            "err_msgs": err_lst,
            "export_level": export_level,
            "export_style": export_style,
            "output_file_name": output_name,
            "output_file_type": file_type,
            "output_file_data": data_encoded,
            "converted_html": converted_html,
            "not_converted_html": not_converted_html,
        }
    else:
        render_data_dct = {
            "request": request,
            "err_msgs": err_lst,
        }

    return templates.TemplateResponse("converter.html", render_data_dct)


@router.get("/equalizer/", include_in_schema=False)
async def equalizer(request: Request):
    return templates.TemplateResponse(
        "equalizer.html", {"request": request, "out_dct": {}}
    )


@router.post("/equalizer/file/", include_in_schema=False)
async def equalize_file(
    request: Request,
    file_obj: UploadFile = File(...),
    match_levels: str = Form(...)
):
    table_info, err_lst = get_table(file_obj, err_lst=[])
    if table_info:
        input_data = InputDictData(data=table_info)
        usr_levels = get_levels(match_levels)
        export_data = await api.equalize_dict(input_data, usr_levels)
        data_encoded = get_url_safe_str(export_data.dict().get("data"))
        file_type = ".xlsx"
        output_name = get_output_name("Equalizer", file_type)
        render_data_dct = {
            "request": request,
            "err_msgs": err_lst,
            "output_file_name": output_name,
            "output_file_type": file_type,
            "output_file_data": data_encoded,
        }
    else:
        render_data_dct = {
            "request": request,
            "err_msgs": err_lst,
        }

    return templates.TemplateResponse("equalizer.html", render_data_dct)


@router.get(
    "/downloads/{data}/{file_type}/{file_name}",
    name="get_download_file",
    include_in_schema=False,
)
def get_download_file(data: str, file_type: str, file_name: str):
    decoded_data = base64.urlsafe_b64decode(data.encode("utf-8"))
    data = json.loads(decoded_data.decode("utf-8"))
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
        return StreamingResponse(excel_io, media_type="application/vnd.ms-excel")


@router.get(f"/", include_in_schema=False)
async def home(request: Request):
    return templates.TemplateResponse(
        "home.html",
        {"request": request, "lynx_version": lynx_version, "api_version": api_version},
    )


@router.get("/levels", include_in_schema=False)
def levels(request: Request):
    return templates.TemplateResponse("levels.html", {"request": request})


@router.get("/nomenclature", include_in_schema=False)
def nomenclature(request: Request):
    return templates.TemplateResponse("nomenclature.html", {"request": request})


@router.get("/static/images/{image_name}", include_in_schema=False)
def get_img(image_name: str):

    return


@router.get("/user_guide", include_in_schema=False)
def user_guide(request: Request):
    return templates.TemplateResponse(
        "user_guide.html",
        {"request": request, "lynx_version": lynx_version, "api_version": api_version},
    )
