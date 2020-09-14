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

import json
import time

import zmq

from lynx.controllers.converter import Converter
from lynx.models.api_models import ConverterExportData
from lynx.utils.toolbox import get_style_level
from lynx.utils.file_handler import (
    table2html,
)
from lynx.utils.job_manager import save_session


def default_worker(worker_id: int):
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.connect("tcp://localhost:5560")

    while True:
        message = socket.recv()
        print(f"Worker #{worker_id} Received request: {message}")
        print(message)
        time.sleep(60)
        try:
            msg = json.loads(message.decode())
        except:
            msg = {"Error": "Error", "Passed": False}
        msg["results"] = f"Finished by worker #{worker_id}"
        msg["time"] = time.time()
        socket.send(json.dumps(msg).encode())


def converter_worker(worker_id: int):
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.connect("tcp://localhost:5560")
    print(f"Worker#{worker_id} started.")
    while True:
        message = socket.recv()
        print(f"Worker #{worker_id} Received Job: {message}")
        try:
            data = json.loads(message.decode()).get("data")
            token = json.loads(message.decode()).get("token", "Temp_token")
            # input_data = InputListData(data=data.get("names", []))
            data_lst = data.get("names", [])
            export_style = data.get("export_style")
            export_level = data.get("export_level")
            style, level = get_style_level(export_style, export_level)
            lynx_converter = Converter(style=style)
            converted_results = ConverterExportData(
                data={"TextInput": lynx_converter.convert_list(data_lst, level=level)}
            )
            # converted_html, not_converted_html = table2html(converted_results)
            response_data = {
                "token": token,
                "err_msgs": [],
                "export_level": export_level,
                "export_style": export_style,
                "converted_results": converted_results.dict(),
                # "converted_html": converted_html,
                # "not_converted_html": not_converted_html,
            }
            save_session(token, response_data)
        except Exception as e:
            response_data = {"Error": e, "Passed": False}

        socket.send(json.dumps(response_data).encode())
