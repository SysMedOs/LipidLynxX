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
from lynx.models.api_models import ConverterExportData, JobType
from lynx.routers.api_linker import link_lipid
from lynx.utils.cfg_reader import app_prefix
from lynx.utils.file_handler import create_converter_output
from lynx.utils.job_manager import save_job
from lynx.utils.temp_file_cleaner import remove_temp_file
from lynx.utils.toolbox import get_style_level, get_url_safe_str


def run_convert_list(token: str, data: dict) -> dict:
    data_lst = data.get("data", [])
    params = data.get("params", {})
    export_style = params.get("style", "ShortHand")
    export_level = params.get("level", "B1")
    file_type = params.get("file_type", "xlsx")
    style, level = get_style_level(export_style, export_level)
    time_tag = time.strftime("%Y%m%d-%H%M%S")
    export_name = f"LipidLynxX-Converter-{time_tag}-{token[:4]}.{file_type}"
    lynx_converter = Converter(style=style)
    converter_results = lynx_converter.convert_list(data_lst, level=level)
    converted_results = ConverterExportData(data={"TextInput": converter_results})
    export_url = f"{app_prefix}/downloads/{export_name}"
    response_data = {
        "token": token,
        "err_msgs": [],
        "export_level": export_level,
        "export_style": export_style,
        "export_name": export_name,
        "export_url": export_url,
        "results": converted_results.dict(),
    }
    save_job(token, response_data)

    return response_data


def run_convert_dict(token: str, data: dict) -> dict:
    data_dct = data.get("data", {})
    params = data.get("params", {})
    export_style = params.get("style", "ShortHand")
    export_level = params.get("level", "B1")
    file_type = params.get("file_type", "xlsx")
    style, level = get_style_level(export_style, export_level)

    results_dct = {}
    for col_name in data_dct:
        print(col_name)
        data_lst = data_dct[col_name]
        lynx_converter = Converter(style=style)
        converter_results = lynx_converter.convert_list(data_lst, level=level)
        results_dct[col_name] = converter_results
    time_tag = time.strftime("%Y%m%d-%H%M%S")
    export_name = f"LipidLynxX-Converter-{time_tag}-{token[:4]}.{file_type}"
    export_url = f"{app_prefix}/downloads/{export_name}"
    converted_results = ConverterExportData(data=results_dct)
    response_data = {
        "token": token,
        "err_msgs": [],
        "export_level": export_level,
        "export_style": export_style,
        "export_name": export_name,
        "export_url": export_url,
        "results": converted_results.dict(),
    }
    save_job(token, response_data)

    return response_data


def run_equalize_list(token: str, data: dict) -> dict:
    pass


def run_link_list(token: str, data: dict) -> dict:

    data_lst = data.get("data", [])
    params = data.get("params", {})
    file_type = params.get("file_type", "xlsx")
    time_tag = time.strftime("%Y%m%d-%H%M%S")
    export_name = f"LipidLynxX-Linker-{time_tag}-{token[:4]}.{file_type}"
    export_url = f"{app_prefix}/downloads/{export_name}"
    all_resources = {}
    export_file_data = {}
    lynx_names = {}
    for lipid_name in data_lst:
        resource_info = link_lipid(lipid_name, export_url=True)
        all_resources[lipid_name] = get_url_safe_str(resource_info)
        export_file_data[lipid_name] = resource_info
        lynx_names[lipid_name] = resource_info.get("lynx_name", "")
    response_data = {
        "token": token,
        "err_msgs": [],
        "export_name": export_name,
        "export_url": export_url,
        "results": {
            "Text input": {
                "all_resources": all_resources,
                "export_file_data": export_file_data,
                "lynx_names": lynx_names,
            }
        },
    }
    save_job(token, response_data)

    return response_data


def init_runner(job: JobType, token: str, data: dict, job_data_type: str):

    if job.job == "convert":
        if job_data_type == "list":
            response_data = run_convert_list(token, data)
        elif job_data_type == "dict":
            print("data")
            print(data)
            response_data = run_convert_dict(token, data)
        else:
            response_data = run_convert_dict(token, data)
    elif job.job == "equalize":
        response_data = run_equalize_list(token, data)
    elif job.job == "link":
        response_data = run_link_list(token, data)
    else:
        response_data = {"token": token, "err_msgs": []}
    remove_temp_file()

    return response_data


def general_worker(worker_id: int, zmq_worker_port: int = 2410):
    context = zmq.Context()
    # Static analysis of Pycharm/PyLint cannot find runtime-defined names e.g. REP
    # See: https://github.com/zeromq/pyzmq/issues/1018
    socket = context.socket(zmq.REP)
    socket.connect(f"tcp://localhost:{zmq_worker_port}")
    print(f"Worker#{worker_id} started @[{zmq_worker_port}].")
    while True:
        message = socket.recv()
        print(f"Worker #{worker_id} @[{zmq_worker_port}] Received Job: {message}")
        try:
            msg_dct = json.loads(message.decode())
            data = msg_dct.get("data")
            token = msg_dct.get("token", "Temp_token")
            job_name = msg_dct.get("job", "").lower()
            job_data_type = msg_dct.get("data_type", "").lower()
            print("job_data_type", job_data_type)
            job = JobType(job=job_name)
            if job.job and token and isinstance(data, dict):
                response_data = init_runner(job, token, data, job_data_type)
            else:
                response_data = {
                    "token": token,
                    "err_msgs": ["Cannot load job information."],
                }
        except Exception as e:
            response_data = {
                "token": "unknown",
                "err_msgs": ["Cannot load job information.", str(e)],
            }

        socket.send(json.dumps(response_data).encode())
