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

import webbrowser
import uvicorn

from lynx.app import app, app_url, app_port


def start_lynx():
    print("Start Browser: ")
    webbrowser.open(f"http://{app_url}:{app_port}", new=1, autoraise=True)
    print("Start LipidLynxX Server: ")
    uvicorn.run(app, host=app_url, port=app_port)


if __name__ == "__main__":
    start_lynx()
