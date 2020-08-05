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

import base64
import json
import os

from fastapi import APIRouter, File, Form, Request, UploadFile, HTTPException
from fastapi.responses import FileResponse
from fastapi.templating import Jinja2Templates

from lynx.models.api_models import (
    FileType,
    StyleType,
    InputListData,
    InputDictData,
    EqualizerExportData,
    LevelsData,
)
from lynx.models.defaults import default_temp_folder
import lynx.routers.api as api
from lynx.utils.cfg_reader import api_version, lynx_version
from lynx.utils.file_handler import (
    create_converter_output,
    create_equalizer_output,
    create_linker_output,
    get_table,
    get_file_type,
    get_output_name,
    table2html,
)
from lynx.utils.frontend_tools import get_converter_response_data, get_linker_response_data
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
    input_data = InputListData(data=names)
    export_style, export_level = get_style_level(export_style, export_level)
    converted_data = await api.convert_list(export_style, input_data, export_level)
    converted_html, not_converted_html = table2html(converted_data)
    response_data = {
        "request": request,
        "err_msgs": [],
        "export_level": export_level,
        "export_style": export_style,
        "converted_html": converted_html,
        "not_converted_html": not_converted_html,

    }
    response_data = get_converter_response_data(converted_data.dict(), file_type, response_data)

    return templates.TemplateResponse("converter.html", response_data)


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
        response_data = {
            "request": request,
            "err_msgs": [],
            "export_level": export_level,
            "export_style": export_style,
            "converted_html": converted_html,
            "not_converted_html": not_converted_html,

        }
        response_data = get_converter_response_data(converted_data.dict(), file_type, response_data)

    else:
        response_data = {
            "request": request,
            "err_msgs": err_lst,
            "output_generated": False,
        }

    return templates.TemplateResponse("converter.html", response_data)


@router.get("/equalizer/", include_in_schema=False)
async def equalizer(request: Request):
    return templates.TemplateResponse(
        "equalizer.html", {"request": request, "out_dct": {}}
    )


@router.post("/equalizer/file/", include_in_schema=False)
async def equalizer_file(
    request: Request, file_obj: UploadFile = File(...), match_levels: str = Form(...)
):
    table_info, err_lst = get_table(file_obj, err_lst=[])
    if table_info:
        input_data = InputDictData(data=table_info)
        usr_levels = get_levels(match_levels)
        if isinstance(usr_levels, str):
            export_data = await api.equalize_single_level(input_data, usr_levels)
        elif isinstance(usr_levels, list):
            levels_data = LevelsData(levels=usr_levels)
            export_data = await api.equalize_multiple_levels(input_data, levels_data)
        else:
            export_data = None
            err_lst.append(f"Invalid levels: {match_levels}")
        if isinstance(export_data, EqualizerExportData):
            # data_encoded = get_url_safe_str(export_data.dict().get("data"))
            file_type = "xlsx"
            output_name = get_output_name("Equalizer", file_type)
            output_path = os.path.join(default_temp_folder, output_name)
            output_info = create_equalizer_output(
                export_data.dict().get("data"), output_name=output_path
            )
            if (
                isinstance(output_info, str)
                and output_info == output_path
                and os.path.isfile(output_path)
            ):
                render_data_dct = {
                    "request": request,
                    "err_msgs": err_lst,
                    "output_file_name": output_name,
                    "output_generated": True,
                }
            else:
                render_data_dct = {
                    "request": request,
                    "err_msgs": err_lst,
                    "output_generated": False,
                }
        else:
            render_data_dct = {
                "request": request,
                "err_msgs": err_lst,
                "output_generated": False,
            }
    else:
        render_data_dct = {
            "request": request,
            "err_msgs": err_lst,
            "output_generated": False,
        }

    return templates.TemplateResponse("equalizer.html", render_data_dct)


@router.get("/linker/", include_in_schema=False)
async def linker(request: Request,):
    return templates.TemplateResponse(
        "linker.html", {"request": request, "all_resources": {}}
    )


@router.post("/linker/list", include_in_schema=False)
async def linker_list(
    request: Request,
    lipid_names: str = Form(...),
    export_url: str = Form(...),
    file_type: FileType = Form(...),
):
    if not lipid_names:
        raise HTTPException(status_code=404)
    if export_url == "include":
        export_url = True
    else:
        export_url = False
    names = lipid_names.splitlines()
    all_resources = {}
    export_file_data = {}
    lynx_names = {}
    for lipid_name in names:
        resource_info = await api.link_str(lipid_name, export_url=True)
        all_resources[lipid_name] = get_url_safe_str(resource_info.get("resource_data"))
        export_file_data[lipid_name] = resource_info
        lynx_names[lipid_name] = resource_info.get("lynx_name", "")

    response_data = {
        "request": request,
        "export_url": export_url,
        "all_resources": all_resources,
        "lynx_names": lynx_names,
    }
    response_data = get_linker_response_data(export_file_data, file_type, response_data)

    return templates.TemplateResponse("linker.html", response_data)


@router.post("/linker/file", include_in_schema=False)
async def linker_file(request: Request,):
    return templates.TemplateResponse(
        "linker.html", {"request": request, "out_dct": {}}
    )


@router.post("/linker/lipid/", include_in_schema=False)
async def linker_lipid(request: Request, lipid_name: str = Form(...)):
    resource_info = await api.link_str(lipid_name, export_url=True)
    resource_info["request"] = request
    return templates.TemplateResponse("resources.html", resource_info)


@router.get("/linker/resources/{data}", include_in_schema=False)
async def view_lipid_resources(request: Request, data: str):
    decoded_data = base64.urlsafe_b64decode(data.encode("utf-8"))
    resource_info = json.loads(decoded_data.decode("utf-8"))
    resource_info["request"] = request
    return templates.TemplateResponse("resources.html", resource_info)


@router.get("/downloads/{file_name}", name="get_download_file", include_in_schema=False)
async def get_download_file(file_name: str):
    file_path = os.path.join(default_temp_folder, file_name)
    if os.path.isfile(file_path):
        response = FileResponse(file_path)
        response.headers["Content-Disposition"] = f"attachment; filename={file_name}"
        return response
    else:
        return f"Failed to generate output file: {file_name}"


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


@router.get("/user_guide", include_in_schema=False)
def user_guide(request: Request):
    return templates.TemplateResponse(
        "guide.html",
        {"request": request, "lynx_version": lynx_version, "api_version": api_version},
    )
