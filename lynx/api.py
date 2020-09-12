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

from fastapi import FastAPI
from fastapi.openapi.utils import get_openapi

from lynx.routers import converter, equalizer, linker  # , parser
from lynx.utils.cfg_reader import api_version, lynx_version

api = FastAPI()
# load APIs from sub-routers
api.include_router(converter.router, prefix="/converter", tags=["converter"])
api.include_router(equalizer.router, prefix="/equalizer", tags=["equalizer"])
api.include_router(linker.router, prefix="/linker", tags=["linker"])
# api.include_router(parser.router, prefix="/parser", tags=["parser"])
# api.include_router(jobs.jobs, prefix="/jobs", tags=["jobs"])


def custom_openapi():
    if api.openapi_schema:
        return api.openapi_schema
    openapi_schema = get_openapi(
        title="LipidLynxX API",
        version=api_version,
        description=f"This is the api (V{api_version}) used in LipidLynxX (V{lynx_version})",
        routes=api.routes,
    )
    openapi_schema["info"]["x-logo"] = {"url": "images/LipidLynxX_icon.png"}
    api.openapi_schema = openapi_schema
    return api.openapi_schema


api.openapi = custom_openapi
