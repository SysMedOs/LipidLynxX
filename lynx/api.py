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

from lynx.routers import api_converter, api_equalizer, api_linker  # , parser
from lynx.utils.cfg_reader import api_version, lynx_version

# create API as a sub-app for lynx.app.py
api = FastAPI(
    debug=True,
    title="LipidLynxX API",
    description=f"This is the api (V{api_version}) used in LipidLynxX (V{lynx_version})",
)
# load APIs from sub-routers
api.include_router(api_converter.router, prefix="/converter", tags=["converter"])
api.include_router(api_equalizer.router, prefix="/equalizer", tags=["equalizer"])
api.include_router(api_linker.router, prefix="/linker", tags=["linker"])
# api.include_router(parser.router, prefix="/parser", tags=["parser"])
# api.include_router(jobs.jobs, prefix="/jobs", tags=["jobs"])

if __name__ == "__main__":
    import uvicorn

    from lynx.mq import start_zmq
    from lynx.utils.cfg_reader import app_cfg_info

    # Start message queue powered by ZeroMQ.
    start_zmq()

    # load app_url and app_port here
    app_url = app_cfg_info.get("app_url", "127.0.0.1")
    app_port = int(app_cfg_info.get("app_port", 1399))

    print("Start LipidLynxX Main Application...")
    uvicorn.run(api, host=app_url, port=app_port)
