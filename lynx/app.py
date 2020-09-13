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
from fastapi.staticfiles import StaticFiles

from lynx.api import api
from lynx.routers.frontend import frontend
from lynx.utils.cfg_reader import app_cfg_info, api_version, lynx_version

# load app_url and app_port here to simplify the LipidLynxX.py file
app_url = app_cfg_info.get("app_url", "127.0.0.1")
app_port = int(app_cfg_info.get("app_port", 1399))

app = FastAPI(debug=True, title="LipidLynxX")

app.mount("/images", StaticFiles(directory="lynx/static/images"), name="images")
# load lynx.api.py as sub-app to provide API service for the frontend
# load frontend from lynx.router.frontend to provide GUI for users
app.include_router(frontend)
# api.py can be started separately to provide API service only
app.mount("/api", api, name="api")


# def custom_openapi():
#     """
#     Settings for the API documentation.
#     """
#     if api.openapi_schema:
#         return api.openapi_schema
#     openapi_schema = get_openapi(
#         title="LipidLynxX API",
#         version=api_version,
#         description=f"This is the api (V{api_version}) used in LipidLynxX (V{lynx_version})",
#         routes=api.routes,
#     )
#     openapi_schema["info"]["x-logo"] = {"url": "images/LipidLynxX_icon.png"}
#     api.openapi_schema = openapi_schema
#     return api.openapi_schema
#
#
# api.openapi = custom_openapi


if __name__ == "__main__":
    import uvicorn

    from lynx.mq import start_zmq

    # Start message queue powered by ZeroMQ.
    start_zmq()

    print("Start LipidLynxX Main Application...")
    uvicorn.run(app, host=app_url, port=app_port)
