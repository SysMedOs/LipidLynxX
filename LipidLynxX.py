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
from lynx.utils.cfg_reader import app_prefix, app_url, app_port


def start_lynx():

    print("Start Browser: ")
    webbrowser.open(
        f"http://{app_url}:{app_port}{app_prefix}", new=1, autoraise=True
    )  # launch default web browser
    # Start message queue powered by ZeroMQ.
    daemon_lynx()
    print("Start LipidLynxX Application...")
    uvicorn.run(app, host=app_url, port=app_port)


if __name__ == "__main__":
    start_lynx()
