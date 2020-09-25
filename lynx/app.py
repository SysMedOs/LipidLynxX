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

# from fastapi.openapi.utils import get_openapi
from fastapi.staticfiles import StaticFiles

from lynx.api import api
from lynx.routers.app_frontend import frontend
from lynx.utils.cfg_reader import app_prefix, app_cfg_info
from lynx.utils.ports import check_port


app = FastAPI(title="LipidLynxX", debug=True)

app.mount(
    f"{app_prefix}/images", StaticFiles(directory="lynx/static/images"), name="images"
)
# load lynx.api.py as sub-app to provide API service for the frontend
# load frontend from lynx.router.frontend to provide GUI for users
app.include_router(frontend, prefix=f"{app_prefix}")
# api.py can be started separately to provide API service only
app.mount(f"{app_prefix}/api", api, name="api")


if __name__ == "__main__":
    import uvicorn

    from lynx.daemon import daemon_lynx

    # Start message queue powered by ZeroMQ.
    # check ports
    checked_zmq_client_port = check_port(
        int(app_cfg_info.get("zmq_client_port", 2409)), task_name="ZMQ client"
    )
    checked_zmq_worker_port = check_port(
        int(app_cfg_info.get("zmq_worker_port", 2410)), task_name="ZMQ worker"
    )
    # run zmq daemon
    daemon_lynx(checked_zmq_client_port, checked_zmq_worker_port)
    daemon_lynx()

    print("Start LipidLynxX Main Application...")
    # check port
    app_url = app_cfg_info.get("app_url", "127.0.0.1")
    app_port = int(app_cfg_info.get("app_port", 1399))
    checked_app_port = check_port(app_port, task_name="LipidLynxX main app")
    if app_port != int(checked_app_port):
        print(f"Port: [{app_port}] in config.ini is already in use.")
        app_port = int(checked_app_port)
        print(f"[INFO] LipidLynxX is now running on port [{checked_app_port}].")
    uvicorn.run(app, host=app_url, port=app_port)
