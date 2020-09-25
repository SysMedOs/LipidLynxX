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
import os
import time

import zmq

from lynx.controllers.converter import Converter
from lynx.models.api_models import ConverterExportData, JobType
from lynx.models.defaults import default_temp_folder
from lynx.utils.cfg_reader import app_prefix
from lynx.utils.file_handler import create_converter_output, table2html
from lynx.utils.job_manager import save_job
from lynx.utils.temp_file_cleaner import remove_temp_file
from lynx.utils.toolbox import get_style_level


# def default_worker(worker_id: int):
#     context = zmq.Context()
#     socket = context.socket(zmq.REP)
#     socket.connect(f"tcp://localhost:{default_zmq_worker_port}")
#
#     while True:
#         message = socket.recv()
#         print(f"Worker #{worker_id} Received request: {message}")
#         print(message)
#         time.sleep(60)
#         try:
#             msg = json.loads(message.decode())
#         except:
#             msg = {"Error": "Error", "Passed": False}
#         msg["results"] = f"Finished by worker #{worker_id}"
#         msg["time"] = time.time()
#         socket.send(json.dumps(msg).encode())


def run_convert_list(token: str, data: dict, zmq_worker_port: int = 2410) -> dict:
    data_lst = data.get("data", [])
    params = data.get("params", {})
    export_style = params.get("style", "ShortHand")
    export_level = params.get("level", "B1")
    file_type = params.get("file_type", "csv")
    style, level = get_style_level(export_style, export_level)
    time_tag = time.strftime("%Y%m%d-%H%M%S")
    export_name = f"LipidLynxX-Converter-{time_tag}-{token[:4]}.{file_type}"
    lynx_converter = Converter(style=style)
    converter_results = lynx_converter.convert_list(data_lst, level=level)
    converted_results = ConverterExportData(data={"TextInput": converter_results})
    output_path = os.path.join(default_temp_folder, export_name)
    export_abs_path = create_converter_output(
        converter_results.dict(), output_name=output_path
    )
    export_url = f"{app_prefix}/downloads/{export_name}"
    # converted_html, not_converted_html = table2html(converted_results)
    response_data = {
        "token": token,
        "err_msgs": [],
        "export_level": export_level,
        "export_style": export_style,
        "export_name": export_name,
        "export_url": export_url,
        "converted_results": converted_results.dict(),
        # "converted_html": converted_html,
        # "not_converted_html": not_converted_html,
    }
    save_job(token, response_data)

    return response_data


def general_worker(worker_id: int, zmq_worker_port: int = 2410):
    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.connect(f"tcp://localhost:{zmq_worker_port}")
    print(f"Worker#{worker_id} started.")
    while True:
        message = socket.recv()
        print(f"Worker #{worker_id} @[{zmq_worker_port}] Received Job: {message}")
        try:
            data = json.loads(message.decode()).get("data")
            token = json.loads(message.decode()).get("token", "Temp_token")
            job_name = json.loads(message.decode()).get("job", "").lower()
            job = JobType(job=job_name)
            # input_data = InputListData(data=data.get("names", []))
            if job.job == "convert":
                response_data = run_convert_list(token, data)
            else:
                response_data = {"token": token, "err_msgs": []}
            remove_temp_file()
        except Exception as e:
            response_data = {
                "token": "unknown",
                "err_msgs": ["Cannot load job information.", str(e)],
            }

        socket.send(json.dumps(response_data).encode())
