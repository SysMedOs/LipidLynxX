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

from fastapi import (
    APIRouter,
    File,
    Form,
    Request,
    UploadFile,
    HTTPException,
)
from fastapi.responses import FileResponse
from fastapi.templating import Jinja2Templates

from lynx.models.api_models import (
    FileType,
    StyleType,
    InputListData,
    InputDictData,
    EqualizerExportData,
    LevelsData,
    JobStatus,
)
from lynx.models.defaults import default_temp_folder, default_template_files
import lynx.api as api
from lynx.routers.api_converter import create_convert_dict_job, create_convert_list_job
from lynx.routers.api_equalizer import create_equalize_job
from lynx.routers.api_linker import create_link_dict_job, create_link_list_job
from lynx.utils.cfg_reader import api_version, lynx_version, app_prefix
from lynx.utils.file_handler import (
    create_equalizer_output,
    get_file_type,
    get_table,
    get_output_name,
)
from lynx.utils.toolbox import get_levels, get_style_level, get_url_safe_str

# upload file size check
# async def valid_content_length(content_length: int = Header(..., lt=10240_000)):
#     return content_length


frontend = APIRouter()
templates = Jinja2Templates(directory="lynx/templates")


# TemplateResponse from jinja2, ignored in API docs page

# get methods
# top-level pages
@frontend.get("/about", include_in_schema=False)
def about(request: Request):
    return templates.TemplateResponse(
        "about.html",
        {"request": request, "lynx_version": lynx_version, "api_version": api_version},
    )


@frontend.get("/converter/", include_in_schema=False)
async def converter(request: Request):
    return templates.TemplateResponse(
        "converter.html", {"request": request, "out_dct": {}}
    )


@frontend.get("/equalizer/", include_in_schema=False)
async def equalizer(request: Request):
    return templates.TemplateResponse(
        "equalizer.html", {"request": request, "out_dct": {}}
    )


@frontend.get(f"/", include_in_schema=False)
async def home(request: Request):
    return templates.TemplateResponse(
        "home.html",
        {"request": request, "lynx_version": lynx_version, "api_version": api_version},
    )


@frontend.get("/levels", include_in_schema=False)
def levels(request: Request):
    return templates.TemplateResponse("levels.html", {"request": request})


@frontend.get("/linker/", include_in_schema=False)
async def linker(request: Request,):
    return templates.TemplateResponse(
        "linker.html", {"request": request, "all_resources": {}}
    )


@frontend.get("/nomenclature", include_in_schema=False)
def nomenclature(request: Request):
    return templates.TemplateResponse("nomenclature.html", {"request": request})


@frontend.get("/user-guide", include_in_schema=False)
def user_guide(request: Request):
    return templates.TemplateResponse(
        "guide.html",
        {"request": request, "lynx_version": lynx_version, "api_version": api_version},
    )


# other functions
@frontend.get(
    "/downloads/{file_name}", name="get_download_file", include_in_schema=False
)
async def get_download_file(file_name: str):
    if file_name in default_template_files:
        file_path = default_template_files.get(file_name)
        path_head, file_name = os.path.split(file_path)
    else:
        file_path = os.path.join(default_temp_folder, file_name)
    print(file_path, file_name)
    if os.path.isfile(file_path):
        response = FileResponse(file_path, filename=file_name)
        # response.headers["Content-Disposition"] = f"attachment; filename={file_name}"
        return response
    else:
        return f"Failed to generate output file: {file_name}"


# post methods
@frontend.post("/converter/text/results/", include_in_schema=False)
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
    export_file_type = get_file_type(file_type)
    job_info = await create_convert_list_job(
        data=input_data,
        style=export_style,
        level=export_level,
        file_type=export_file_type,
    )  # type: JobStatus
    response_data = job_info.dict()
    response_data["request"] = request
    # rel link to get_download
    response_data["export_level"] = export_level
    response_data["export_style"] = export_style
    response_data["err_msgs"] = []

    return templates.TemplateResponse("converter_results.html", response_data)


@frontend.post("/converter/file/results/", include_in_schema=False)
async def converter_file(
    request: Request,
    file_obj: UploadFile = File(...),
    export_level: str = Form(...),
    export_style: StyleType = Form(...),
    file_type: FileType = Form(...),
):
    table_info, err_lst = get_table(file_obj, err_lst=[])
    export_file_type = get_file_type(file_type)
    export_style, export_level = get_style_level(export_style, export_level)
    if table_info:
        input_data = InputDictData(data=table_info)
        job_info = await create_convert_dict_job(
            data=input_data,
            style=export_style,
            level=export_level,
            file_type=export_file_type,
        )  # type: JobStatus
        response_data = job_info.dict()
        response_data["err_msgs"] = []
    else:
        response_data = {"err_msgs": [f"Failed to read file: {file_obj.filename}"]}
    response_data["request"] = request
    response_data["export_level"] = export_level
    response_data["export_style"] = export_style

    return templates.TemplateResponse("converter_results.html", response_data)


@frontend.post("/equalizer/file/", include_in_schema=False)
async def equalizer_file(
    request: Request, file_obj: UploadFile = File(...), match_levels: str = Form(...)
):
    table_info, err_lst = get_table(file_obj, err_lst=[])

    if table_info:
        input_data = InputDictData(data=table_info)
        usr_levels = get_levels(match_levels)
        if isinstance(usr_levels, str):
            usr_levels = [usr_levels]
        job_info = await create_equalize_job(
            data=input_data, levels=usr_levels, file_type="xlsx"
        )  # type: JobStatus
        response_data = job_info.dict()
        response_data["err_msgs"] = []
    else:
        response_data = {"err_msgs": [f"Failed to read file: {file_obj.filename}"]}
    response_data["request"] = request

    return templates.TemplateResponse("equalizer_results.html", response_data)


@frontend.post("/linker/text", include_in_schema=False)
async def linker_text(
    request: Request, lipid_names: str = Form(...), file_type: FileType = Form(...),
):
    if not lipid_names:
        raise HTTPException(status_code=404)

    names = lipid_names.splitlines()
    export_file_type = get_file_type(file_type)
    input_data = InputListData(data=names)
    job_info = await create_link_list_job(
        data=input_data, file_type=export_file_type,
    )  # type: JobStatus
    response_data = job_info.dict()
    response_data["request"] = request
    # rel link to get_download
    response_data["err_msgs"] = []

    return templates.TemplateResponse("linker_results_summary.html", response_data)


@frontend.post("/linker/file", include_in_schema=False)
async def linker_file(
    request: Request, file_obj: UploadFile = File(...), file_type: FileType = Form(...),
):
    table_info, err_lst = get_table(file_obj, err_lst=[])
    export_file_type = get_file_type(file_type)
    print(table_info)
    if table_info:
        input_data = InputDictData(data=table_info)
        job_info = await create_link_dict_job(
            data=input_data, file_type=export_file_type,
        )  # type: JobStatus
        response_data = job_info.dict()
        response_data["request"] = request
        # rel link to get_download
        response_data["err_msgs"] = []
        # resource_info = await api.link_dict(table_info, export_url=True)
        # for k in resource_info:
        #     if len(k) > 16:
        #         display_k = f"{k[:15]}~{k[-1]}"
        #     else:
        #         display_k = k
        #     col_resource_info = resource_info.get(k, [])
        #     lynx_names = {}
        #     all_resources = {}
        #     for lipid_name in col_resource_info:
        #         lipid_name_info = col_resource_info.get(lipid_name, {})
        #         all_resources[lipid_name] = get_url_safe_str(lipid_name_info)
        #         lynx_names[lipid_name] = lipid_name_info.get("lynx_name")
        #     sum_all_resources[display_k] = {
        #         "all_resources": all_resources,
        #         "export_file_data": col_resource_info,
        #         "lynx_names": lynx_names,
        #     }
    else:
        response_data = {"err_msgs": [f"Failed to read file: {file_obj.filename}"]}
    response_data["request"] = request
    return templates.TemplateResponse("linker_results_summary.html", response_data)


# direct fast link of one lipid on the home page
@frontend.post("/linker/results/details", include_in_schema=False)
async def linker_lipid(
    request: Request,
    lipid_name: str = Form(...),
    resource_data: str = Form(...),
    export_name: str = Form(...),
    export_url: str = Form(...),
):
    if resource_data == "NO_RESOURCE_DATA":
        resource_info = await api.api_linker.get_link_lipid(lipid_name, export_url=True)
    else:
        decoded_data = base64.urlsafe_b64decode(resource_data.encode("utf-8"))
        resource_info = json.loads(decoded_data.decode("utf-8"))
    resource_info["request"] = request
    if export_name.lower().endswith("xlsx") and export_url.lower().endswith("xlsx"):
        resource_info["export_name"] = export_name
        resource_info["export_url"] = export_url
    return templates.TemplateResponse("linker_results_details.html", resource_info)


# @frontend.post("/linker/results/summary", include_in_schema=False)
# async def view_resource_summary(
#     request: Request,
#     display_data: str = Form(...),
#     lynx_names: str = Form(...),
#     file_name: str = Form(...),
#     export_url: bool = Form(...),
# ):
#     decoded_data = base64.urlsafe_b64decode(display_data.encode("utf-8"))
#     decoded_lynx_names = base64.urlsafe_b64decode(lynx_names.encode("utf-8"))
#     resource_info = {
#         "data": json.loads(decoded_data.decode("utf-8")),
#         "lynx_names": json.loads(decoded_lynx_names.decode("utf-8")),
#         "file_name": file_name,
#         "export_url": export_url,
#         "request": request,
#     }
#     return templates.TemplateResponse("linker_results_summary.html", resource_info)
