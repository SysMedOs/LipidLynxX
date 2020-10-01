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

import os

from lynx.models.api_models import ConverterExportData, FileType
from lynx.models.defaults import default_temp_folder
from lynx.utils.file_handler import (
    create_converter_output,
    create_equalizer_output,
    create_linker_output,
    get_file_type,
    get_output_name,
)
from lynx.utils.toolbox import get_url_safe_str


# def get_response_data(
#     output_info: str, output_path: str, output_name: str, response_data: dict
# ):
#     if (
#         isinstance(output_info, str)
#         and output_info == output_path
#         and os.path.isfile(output_path)
#     ):
#         response_data["output_generated"] = True
#         response_data["output_file_name"] = output_name
#     else:
#         response_data["output_generated"] = False
#         response_data["output_file_name"] = ""
#         response_data["err_msgs"] = response_data["err_msgs"] + [
#             "Failed to generate output."
#         ]
#     return response_data


# def get_converter_response_data(data: dict, file_type: FileType, response_data: dict):
#     file_type = get_file_type(file_type)
#     output_name = get_output_name("Converter", file_type)
#     output_path = os.path.join(default_temp_folder, output_name)
#     output_info = create_converter_output(
#         data.get("data"), output_name=output_path, file_type=file_type
#     )
#     response_data = get_response_data(
#         output_info, output_path, output_name, response_data
#     )
#
#     return response_data
#
#
# def get_linker_response_data(response_data: dict, file_type: FileType = FileType.xlsx):
#     export_data = response_data.get("data", {})
#     if export_data:
#         file_type = get_file_type(file_type)
#         output_name = get_output_name("Linker", file_type)
#         output_path = os.path.join(default_temp_folder, output_name)
#         export_url = response_data.get("export_url", True)
#         output_info = create_linker_output(
#             response_data.get("data", {}),
#             output_name=output_path,
#             file_type=file_type,
#             export_url=export_url,
#         )
#         data = {}
#         for col in export_data:
#             col_data = export_data.get(col, {})
#             all_resources = col_data.get("all_resources", {})
#             lynx_names = col_data.get("lynx_names", {})
#             data[col] = {
#                 "display_data": get_url_safe_str(all_resources),
#                 "lynx_names": get_url_safe_str(lynx_names),
#             }
#         response_data["data"] = data
#         response_data = get_response_data(
#             output_info, output_path, output_name, response_data
#         )
#
#     return response_data
