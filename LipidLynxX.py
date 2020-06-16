# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
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
