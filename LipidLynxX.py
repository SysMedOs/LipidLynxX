# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
import multiprocessing
import webbrowser

import uvicorn

import lynx
from lynx.models.defaults import cfg_info_dct
from lynx.controllers.rest import apilynx


def api_server():
    print("Start API service: ", multiprocessing.current_process().name)
    api_url = cfg_info_dct.get("api_url", "http://127.0.0.1")
    api_port = int(cfg_info_dct.get("api_port", 1399))
    print(f"API {api_url} {api_port}")
    uvicorn.run(apilynx.api_app, host="localhost", port=api_port)


def website():
    print("Start GUI: ", multiprocessing.current_process().name)
    base_url = cfg_info_dct.get("base_url", "http://localhost")
    base_port = int(cfg_info_dct.get("base_port", 1451))
    # webbrowser.open(f"{base_url}:{base_port}/lynx", new=1, autoraise=True)
    lynx.app.run(debug=True, host=base_url, port=base_port)


def browser():
    print("Start Browser: ", multiprocessing.current_process().name)
    base_url = cfg_info_dct.get("base_url", "127.0.0.1")
    base_port = int(cfg_info_dct.get("base_port", 1451))
    webbrowser.open(f"http://{base_url}:{base_port}/lynx", new=1, autoraise=True)


if __name__ == '__main__':
    api = multiprocessing.Process(name='Lynx API', target=api_server)
    web = multiprocessing.Process(name='Lynx GUI', target=website)
    ui = multiprocessing.Process(name='Lynx Browser', target=browser)

    api.start()
    web.start()
    ui.start()
