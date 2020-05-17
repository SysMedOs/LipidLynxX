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

from io import BytesIO

from fastapi import FastAPI, Request, File, Form, UploadFile
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import pandas as pd
import uvicorn

api_app = FastAPI(debug=True)
api_app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")

api_version = "0.2"
lynx_version = "0.5.16"


@api_app.get(f"/")
async def home(request: Request):
    return templates.TemplateResponse("home.html",
                                      {"request": request, "lynx_version": lynx_version, "api_version": api_version})


@api_app.get("/converter/")
async def converter(request: Request):
    return templates.TemplateResponse(
        "fast_converter.html",
        {"request": request})


@api_app.post("/converter/text/")
async def convert_text(request: Request, lipid_names: str = Form(...), export_level: str = Form(...),
                       export_style: str = Form(...),
                       file_type: str = Form(...)):
    return templates.TemplateResponse(
        "fast_converter.html",
        {"request": request, "lipid_names": lipid_names, "export_level": export_level,
         "export_style": export_style, "file_type": file_type}
    )


@api_app.post("/converter/file/")
async def convert_file(request: Request, file_obj: UploadFile = File(...), export_level: str = Form(...),
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
    uvicorn.run(api_app, host="127.0.0.1", port=1399)
