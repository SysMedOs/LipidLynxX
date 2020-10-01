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

import webbrowser
import uvicorn

from lynx.app import app
from lynx.daemon import daemon_lynx
from lynx.utils.cfg_reader import app_prefix, app_cfg_info
from lynx.utils.ports import check_port


def start_lynx(app_url: str, app_port: int):

    print("Start Browser: ")
    app_link = f"http://{app_url}:{app_port}{app_prefix}"
    webbrowser.open(app_link, new=1, autoraise=True)  # launch default web browser
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
    print("Start LipidLynxX Application...")
    uvicorn.run(app, host=app_url, port=app_port)


if __name__ == "__main__":
    # check port
    checked_app_url = app_cfg_info.get("app_url", "127.0.0.1")
    app_port = int(app_cfg_info.get("app_port", 1399))
    checked_app_port = check_port(app_port, task_name="LipidLynxX main app")
    if app_port != checked_app_port:
        print(f"Port: [{app_port}] in config.ini is already in use.")
        print(f"[INFO] LipidLynxX is now running on port [{checked_app_port}].")
    start_lynx(checked_app_url, checked_app_port)
