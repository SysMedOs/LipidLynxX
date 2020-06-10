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

from fastapi import APIRouter
from fastapi import Request, File, Form, UploadFile
from fastapi.templating import Jinja2Templates
import pandas as pd

from lynx.config import api_version, lynx_version
from lynx.models.api_models import StyleName, InputListData
import lynx.routers.api as api


router = APIRouter()
templates = Jinja2Templates(directory="lynx/templates")


# TemplateResponse from jinja2, ignored in API docs page
@router.get(f"/", include_in_schema=False)
async def home(request: Request):
    return templates.TemplateResponse("home.html",
                                      {"request": request, "lynx_version": lynx_version, "api_version": api_version})


@router.get(f"/api_docs", include_in_schema=False)
async def docs(request: Request):
    return templates.TemplateResponse("docs.html",
                                      {"request": request})


@router.get("/converter/", include_in_schema=False)
async def converter(request: Request):
    return templates.TemplateResponse(
        "fast_converter.html",
        {"request": request})


@router.post("/converter/text/", include_in_schema=False)
async def converter_text(request: Request, lipid_names: str = Form(...),
                         export_level: str = Form(...),
                         export_style: StyleName = Form(...),
                         file_type: str = Form(...)):
    names = lipid_names.splitlines()
    input_data = InputListData(lipid_names=names)
    converted_data = api.convert_list(export_style, input_data)
    print(converted_data)
    return templates.TemplateResponse(
        "fast_converter.html",
        {"request": request, "lipid_names": lipid_names, "export_level": export_level,
         "export_style": export_style, "file_type": file_type}
    )


@router.post("/converter/file/", include_in_schema=False)
async def converter_file(request: Request, file_obj: UploadFile = File(...), export_level: str = Form(...),
                         export_style: str = Form(...),
                         file_type: str = Form(...)):
    df = pd.read_csv(file_obj.file)
    print(df)
    file_name = "Uploaded"
    return templates.TemplateResponse(
        "fast_converter.html",
        {"request": request, "file_obj": file_obj, "export_level": export_level,
         "export_style": export_style, "file_type": file_type, "file_name": file_name}
    )
