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

import asyncio
import json
import os
import time

import zmq

from lynx.controllers.converter import Converter
from lynx.models.api_models import ConverterExportData, JobType
from lynx.models.defaults import default_temp_folder
from lynx.routers.api_linker import link_one_lipid
from lynx.utils.cfg_reader import app_prefix
from lynx.utils.file_handler import (
    create_converter_output,
    create_equalizer_output,
    create_linker_output,
)
from lynx.utils.job_manager import save_job
from lynx.utils.temp_file_cleaner import remove_temp_file
from lynx.utils.toolbox import get_style_level, get_url_safe_str


def get_export_file_info(token: str, file_type: str, client_name: str):
    time_tag = time.strftime("%Y%m%d-%H%M%S")
    export_name = f"LipidLynxX-{client_name}-{time_tag}-{token[:4]}.{file_type}"
    export_path = os.path.join(default_temp_folder, export_name)
    export_url = f"{app_prefix}/downloads/{export_name}"

    return export_name, export_path, export_url


def convert_list(data: list, style: str, level: str) -> ConverterExportData:
    lynx_converter = Converter(style=style)
    converter_results = lynx_converter.convert_list(data, level=level)
    converted_results = ConverterExportData(data={"TextInput": converter_results})

    return converted_results


def convert_dict(data: dict, style: str, level: str) -> ConverterExportData:
    results_dct = {}
    for col_name in data:
        print(col_name)
        data_lst = data[col_name]
        lynx_converter = Converter(style=style)
        converter_results = lynx_converter.convert_list(data_lst, level=level)
        results_dct[col_name] = converter_results
    converted_results = ConverterExportData(data=results_dct)

    return converted_results


def run_converter(token: str, data: dict, job_data_type: str) -> dict:
    params = data.get("params", {})
    export_style = params.get("style", "ShortHand")
    export_level = params.get("level", "B1")
    file_type = params.get("file_type", "xlsx")
    style, level = get_style_level(export_style, export_level)
    export_name, export_path, export_url = get_export_file_info(
        token, file_type, client_name="Converter"
    )
    if job_data_type == "list":
        data_to_run = data.get("data", [])
        converted_results = convert_list(data_to_run, style=style, level=level)
    elif job_data_type == "dict":
        data_to_run = data.get("data", {})
        converted_results = convert_dict(data_to_run, style=style, level=level)
    else:
        data_to_run = data.get("data")
        if isinstance(data_to_run, list):
            converted_results = convert_list(data_to_run, style=style, level=level)
        elif isinstance(data_to_run, dict):
            converted_results = convert_dict(data_to_run, style=style, level=level)
        else:
            converted_results = ConverterExportData(data={})
    export_abs_path = create_converter_output(
        converted_results.dict().get("data"), output_name=export_path
    )
    response_data = {
        "token": token,
        "err_msgs": [],
        "export_level": export_level,
        "export_style": export_style,
        "export_name": export_name,
        "export_url": export_url,
        "results": converted_results.dict(),
    }
    if os.path.isfile(export_abs_path):
        pass
    else:
        response_data["err_msgs"] = [f"Cannot create output {export_name}"]
    save_job(token, response_data)

    return response_data


def run_equalizer(token: str, data: dict, job_data_type: str) -> dict:
    data_lst = data.get("data", [])
    params = data.get("params", {})
    file_type = params.get("file_type", "xlsx")
    export_name, export_path, export_url = get_export_file_info(
        token, file_type, client_name="Linker"
    )


async def link_list(data: list, export_path: str, file_type: str):

    all_resources = {}
    export_file_data = {}
    lynx_names = {}
    for lipid_name in data:
        resource_info = await link_one_lipid(lipid_name, export_url=True)
        all_resources[lipid_name] = get_url_safe_str(resource_info)
        export_file_data[lipid_name] = resource_info
        lynx_names[lipid_name] = resource_info.get("lynx_name", "")

    export_data = {
        "text_input": {
            "all_resources": all_resources,
            "export_file_data": export_file_data,
            "lynx_names": lynx_names,
        }
    }
    output_info = create_linker_output(
        export_data, output_name=export_path, file_type=file_type, export_url=True,
    )
    data = {}
    for col in export_data:
        col_data = export_data.get(col, {})
        all_resources = col_data.get("all_resources", {})
        lynx_names = col_data.get("lynx_names", {})
        data[col] = {
            "display_data": get_url_safe_str(all_resources),
            "lynx_names": get_url_safe_str(lynx_names),
        }

    return export_data


async def link_dict(data: dict, export_path: str, file_type: str):

    export_data = {}
    for col in data:
        lipid_col = data[col]
        linked_col_info = {}
        all_resources = {}
        export_file_data = {}
        lynx_names = {}
        for lipid_name in lipid_col:
            resource_info = await link_one_lipid(lipid_name, export_url=True)
            all_resources[lipid_name] = get_url_safe_str(resource_info)
            export_file_data[lipid_name] = resource_info
            lynx_names[lipid_name] = resource_info.get("lynx_name", "")
        export_data[col] = {
            "all_resources": all_resources,
            "export_file_data": export_file_data,
            "lynx_names": lynx_names,
        }

    output_info = create_linker_output(
        export_data, output_name=export_path, file_type=file_type, export_url=True,
    )
    data = {}
    for col in export_data:
        col_data = export_data.get(col, {})
        all_resources = col_data.get("all_resources", {})
        lynx_names = col_data.get("lynx_names", {})
        data[col] = {
            "display_data": get_url_safe_str(all_resources),
            "lynx_names": get_url_safe_str(lynx_names),
        }

    return export_data


def run_linker(token: str, data: dict, job_data_type: str):
    params = data.get("params", {})
    file_type = params.get("file_type", "xlsx")
    export_name, export_path, export_url = get_export_file_info(
        token, file_type, client_name="Linker"
    )
    if job_data_type == "list":
        data_to_run = data.get("data", [])
        loop = asyncio.get_event_loop()
        export_data = loop.run_until_complete(
            link_list(data_to_run, export_path=export_path, file_type=file_type)
        )
    elif job_data_type == "dict":
        data_to_run = data.get("data", {})
        loop = asyncio.get_event_loop()
        export_data = loop.run_until_complete(
            link_dict(data_to_run, export_path=export_path, file_type=file_type)
        )
    else:
        data_to_run = data.get("data")
        if isinstance(data_to_run, list):
            loop = asyncio.get_event_loop()
            export_data = loop.run_until_complete(
                link_list(data_to_run, export_path=export_path, file_type=file_type)
            )
        elif isinstance(data_to_run, dict):
            loop = asyncio.get_event_loop()
            export_data = loop.run_until_complete(
                link_dict(data_to_run, export_path=export_path, file_type=file_type)
            )
        else:
            export_data = {}

    response_data = {
        "token": token,
        "err_msgs": [],
        "export_name": export_name,
        "export_url": export_url,
        "results": export_data,
    }

    if os.path.isfile(export_path):
        pass
    else:
        response_data["err_msgs"] = [f"Cannot create output {export_name}"]
    save_job(token, response_data)

    return response_data


def init_runner(job: JobType, token: str, data: dict, job_data_type: str):

    if job.job == "convert":
        print("run_converter")
        response_data = run_converter(token, data, job_data_type)
    elif job.job == "equalize":
        response_data = run_equalizer(token, data, job_data_type)
    elif job.job == "link":
        response_data = run_linker(token, data, job_data_type)
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

        # msg_dct = json.loads(message.decode())
        # data = msg_dct.get("data")
        # token = msg_dct.get("token", "Temp_token")
        # job_name = msg_dct.get("job", "").lower()
        # job_data_type = msg_dct.get("data_type", "").lower()
        # job = JobType(job=job_name)
        # if job.job and token and isinstance(data, dict):
        #     response_data = init_runner(job, token, data, job_data_type)
        # else:
        #     response_data = {
        #         "token": token,
        #         "err_msgs": ["Cannot load job information."],
        #     }

        try:
            msg_dct = json.loads(message.decode())
            data = msg_dct.get("data")
            token = msg_dct.get("token", "Temp_token")
            job_name = msg_dct.get("job", "").lower()
            job_data_type = msg_dct.get("data_type", "").lower()
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
