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

from fastapi import FastAPI
from fastapi.openapi.utils import get_openapi
from fastapi.staticfiles import StaticFiles

from lynx.config import api_version, lynx_version
from lynx.models.defaults import cfg_info_dct
from lynx.routers import api, frontend

app_url = cfg_info_dct.get("app_url", "127.0.0.1")
app_port = int(cfg_info_dct.get("app_port", 1399))

app = FastAPI(debug=True)


app.include_router(api.router, prefix="/api", tags=["api"])
app.include_router(frontend.router)
app.mount("/images", StaticFiles(directory="lynx/static/images"), name="images")


def custom_openapi(openapi_prefix: str):
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="LipidLynxX API",
        version=api_version,
        description=f"This is the api (V{api_version}) used in LipidLynxX (V{lynx_version})",
        routes=app.routes,
        openapi_prefix=openapi_prefix,
    )
    openapi_schema["info"]["x-logo"] = {"url": "static/images/LipidLynxX_Logo.png"}
    app.openapi_schema = openapi_schema
    return app.openapi_schema


app.openapi = custom_openapi

if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host=app_url, port=app_port)
