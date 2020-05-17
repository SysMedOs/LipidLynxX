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

from enum import Enum
from io import BytesIO
from typing import List

from fastapi import Body, FastAPI, Request, File, Form, UploadFile
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel
import pandas as pd
import uvicorn

from lynx.controllers.converter import Converter
from lynx.models.defaults import cfg_info_dct

api_version = "0.2"
lynx_version = "0.5.16"

api_app = FastAPI(debug=True)
api_app.mount("/static", StaticFiles(directory="lynx/static"), name="static")
templates = Jinja2Templates(directory="lynx/templates")


api_url = cfg_info_dct.get("api_url", "127.0.0.1")
api_port = int(cfg_info_dct.get("api_port", 1399))


# Define objects
class FileType(str, Enum):
    xlsx = "xlsx"
    csv = "csv"


class StyleName(str, Enum):
    lipidlynxx = "LipidLynxX"
    comp_db = "COMP_DB"


class InputListData(BaseModel):
    lipid_names: List[str]

    class Config:
        schema_extra = {
            "example": {
                "lipid_names": ["DHA", "PLPC", "PI 16:0-20:4", "UNKOWN"]
            }
        }


class ExportListData(BaseModel):
    input: List[str]
    output: List[str]
    converted: List[List[str]]
    skipped: List[str]

    class Config:
        schema_extra = {
            "example": {
                "input": [
                    "DHA",
                    "PLPC",
                    "PI 16:0-20:4",
                    "UNKOWN"
                ],
                "output": [
                    "FA22:6",
                    "PC(16:0/18:2)",
                    "PI(16:0_20:4)"
                ],
                "converted": [
                    [
                        "DHA",
                        "FA22:6"
                    ],
                    [
                        "PLPC",
                        "PC(16:0/18:2)"
                    ],
                    [
                        "PI 16:0-20:4",
                        "PI(16:0_20:4)"
                    ]
                ],
                "skipped": ["UNKOWN"]
            }
        }


# Pure APIs
@api_app.post("/convert/list/", response_model=ExportListData)
def convert_list(style: StyleName, data: InputListData):
    """
    Convert a list of lipid names into supported abbreviations
    """
    lynx_converter = Converter(style=style)
    print(data)
    x = lynx_converter.convert_list(data.lipid_names)
    return x


# TemplateResponse from jinja2, ignored in API docs page
@api_app.get(f"/", include_in_schema=False)
async def home(request: Request):
    return templates.TemplateResponse("home.html",
                                      {"request": request, "lynx_version": lynx_version, "api_version": api_version})


@api_app.get(f"/api_docs", include_in_schema=False)
async def docs(request: Request):
    return templates.TemplateResponse("docs.html",
                                      {"request": request})


@api_app.get("/converter/", include_in_schema=False)
async def converter(request: Request):
    return templates.TemplateResponse(
        "fast_converter.html",
        {"request": request})


@api_app.post("/converter/text/", include_in_schema=False)
async def converter_text(request: Request, lipid_names: str = Form(...),
                         export_level: str = Form(...),
                         export_style: StyleName = Form(...),
                         file_type: str = Form(...)):
    names = lipid_names.splitlines()
    l_data = ExportListData(lipid_names=names)
    nl = convert_list(export_style, l_data)
    print(nl)
    return templates.TemplateResponse(
        "fast_converter.html",
        {"request": request, "lipid_names": lipid_names, "export_level": export_level,
         "export_style": export_style, "file_type": file_type}
    )


@api_app.post("/converter/file/", include_in_schema=False)
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


if __name__ == '__main__':
    uvicorn.run(api_app, host=api_url, port=api_port)
