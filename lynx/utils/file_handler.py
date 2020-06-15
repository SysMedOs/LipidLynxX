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

from datetime import datetime
import json
import os
from io import BytesIO
from typing import List, Union

from fastapi import File
import pandas as pd

from lynx.models.api_models import ConverterExportData, ConvertedListData, FileType
import lynx.utils
from lynx.utils.log import logger
from lynx.utils.toolbox import keep_string_only


def get_abs_path(file_path: str) -> str:
    """
    check and get absolute file path from given file
    Args:
        file_path: The relative path to a file

    Returns:
        abs_path: the absolute path of input file

    """

    abs_path = ""

    # Force change to LipidLynxX folder for macOS
    cwd = os.getcwd()
    lynx_utils_dir = os.path.dirname(lynx.utils.__file__)
    lynx_dir = os.path.dirname(lynx_utils_dir)
    project_dit = os.path.dirname(lynx_dir)
    if cwd != project_dit:
        os.chdir(project_dit)
        # print("Change to LipidLynxX project folder", os.getcwd())
    else:
        # print("Working on folder", os.getcwd())
        pass

    if os.path.isdir(file_path):
        abs_path = os.path.abspath(file_path)

    elif os.path.isfile(file_path):
        abs_path = os.path.abspath(file_path)
    else:
        in_file_lst = [
            r"{path}".format(path=file_path),
            r"./{path}".format(path=file_path),
            r"../{path}".format(path=file_path),
            r"../../{path}".format(path=file_path),
            r"../../../{path}".format(path=file_path),
        ]
        for f in in_file_lst:
            if os.path.isdir(f):
                abs_path = os.path.abspath(f)
                break
            elif os.path.isfile(f):
                abs_path = os.path.abspath(f)
                break

    if not abs_path:
        raise FileNotFoundError(
            f"Can not find file: {file_path} from path {os.getcwd()}"
        )

    return abs_path


def get_json(file: str) -> dict:
    file = get_abs_path(file)
    if file.lower().endswith(".json"):
        with open(file) as file_obj:
            js_obj = json.load(file_obj)
            return js_obj
    else:
        raise IOError(f"Input file: {file} is not json file")


def load_folder(folder: str, file_type: str = "") -> List[str]:
    """
     Load all files under given folder, optional with selected file suffix
     Args:
         folder: path of the folder.
         file_type: type of the file, default value is "" for no file type filter

     Returns:
         file_abs_path_lst: the list of files under given folder in absolute path

     """
    abs_path = get_abs_path(folder)
    file_lst = os.listdir(abs_path)
    file_abs_path_lst = [os.path.join(abs_path, x) for x in file_lst]
    if file_type:
        file_abs_path_lst = [
            f for f in file_abs_path_lst if f.lower().endswith(file_type.lower())
        ]
    file_abs_path_lst = [abs_f for abs_f in file_abs_path_lst if os.path.isfile(abs_f)]
    logger.debug(
        f"Fund {file_type} files:\nunder folder: {folder}\nfiles:\n {file_abs_path_lst}"
    )

    return file_abs_path_lst


def create_converter_output(
    data: dict, output_name: str = None, file_type: str = ".xlsx"
) -> Union[BytesIO, str]:
    file_info = None
    converted_df = pd.DataFrame()
    not_converted_df = pd.DataFrame()
    if data:
        not_converted_dct = {}
        df_lst = []
        for k in data:
            if isinstance(data[k], dict):
                k_pairs = data[k].get("converted", [])
                k_not_converted = data[k].get("skipped", [])
                if k_pairs and isinstance(k, str):
                    df_lst.append(pd.DataFrame(k_pairs, columns=[k, f"{k}_converted"]))

                if k_not_converted:
                    not_converted_dct[f"{k}_skipped"] = k_not_converted
            elif isinstance(data[k], ConvertedListData):
                k_pairs = data[k].converted
                if k_pairs and isinstance(k, str):
                    df_lst.append(pd.DataFrame(k_pairs, columns=[k, f"{k}_converted"]))
                k_not_converted = data[k].skipped
                if k_not_converted:
                    not_converted_dct[f"{k}_skipped"] = k_not_converted
            elif isinstance(data[k], list) and k == "converted":
                k_pairs = data.get("converted", [])
                if k_pairs:
                    df_lst.append(
                        pd.DataFrame(k_pairs, columns=["input", f"converted"])
                    )
            elif isinstance(data[k], list) and k == "skipped":
                k_not_converted = data.get("skipped", [])
                if k_not_converted:
                    not_converted_dct[f"skipped"] = k_not_converted

        if df_lst:
            converted_df = pd.concat(df_lst, axis=1)

        if not_converted_dct:
            not_converted_df = pd.DataFrame.from_dict(
                not_converted_dct, orient="index"
            ).T
        if not converted_df.empty:
            if output_name and isinstance(output_name, str):
                try:
                    if file_type.lower().endswith("csv"):
                        converted_df.to_csv(output_name)
                    else:
                        converted_df.to_excel(
                            output_name, sheet_name="converted", index=False
                        )
                    file_info = get_abs_path(output_name)
                except IOError:
                    file_info = (
                        f"[IO error] Cannot create file: {output_name} as output."
                    )
            else:
                file_info = BytesIO()
                if file_type.lower().endswith("csv"):
                    file_info.write(converted_df.to_csv().encode("utf-8"))

                else:
                    output_writer = pd.ExcelWriter(
                        file_info, engine="openpyxl"
                    )  # write to BytesIO instead of file path
                    converted_df.to_excel(
                        output_writer, sheet_name="converted", index=False
                    )
                    if not not_converted_df.empty:
                        not_converted_df.to_excel(
                            output_writer, sheet_name="skipped", index=False
                        )
                    output_writer.save()
                file_info.seek(0)

    return file_info


def create_equalizer_output(
    data: dict, output_name: str = None
) -> Union[BytesIO, str]:
    file_info = None
    is_file_name = False
    if data:
        if output_name and isinstance(output_name, str):
            try:
                table_writer = pd.ExcelWriter(output_name, engine="openpyxl")
                is_file_name = True
                file_info = output_name
            except IOError:
                file_info = f"[IO error] Cannot create file: {output_name} as output."
                return file_info
        else:
            file_info = BytesIO()
            table_writer = pd.ExcelWriter(
                file_info, engine="openpyxl"
            )  # write to BytesIO instead of file path
        equalized_data = data.get("equalized")
        for lv in equalized_data:
            lv_data = equalized_data.get(lv, {})
            for k in lv_data:
                if k.lower().startswith("match"):
                    matched_dct = lv_data[k]
                    if matched_dct:
                        out_matched_df = pd.DataFrame.from_dict(
                            matched_dct, orient="index"
                        )
                        out_matched_df.index.names = [f"ID@Lv_{lv}"]
                        out_matched_df.sort_index().to_excel(
                            table_writer, sheet_name=f"Lv_{lv}_matched"
                        )
                    else:
                        pass
                elif k.lower().startswith("unmatched"):
                    unmatched_dct = lv_data[k]
                    if unmatched_dct:
                        out_unmatched_df = pd.DataFrame.from_dict(
                            unmatched_dct, orient="index"
                        )
                        out_unmatched_df.index.names = [f"ID@Lv_{lv}"]
                        out_unmatched_df.sort_index().to_excel(
                            table_writer, sheet_name=f"Lv_{lv}_unmatched"
                        )
                    else:
                        pass
                else:
                    pass
        skipped_dct = data.get("skipped")
        if skipped_dct:
            pd.DataFrame.from_dict(
                skipped_dct, orient="index"
            ).T.sort_index().to_excel(
                table_writer, sheet_name="skipped", index=False
            )
        else:
            pass

        table_writer.save()
        if not is_file_name:
            file_info.seek(0)
    else:
        file_info = None
    return file_info


def get_table(uploaded_file: File, err_lst: list) -> (dict, list):
    usr_file_name = uploaded_file.filename
    if usr_file_name.lower().endswith("xlsx"):
        table_dct = pd.read_excel(uploaded_file.file).to_dict(orient="list")
    elif usr_file_name.lower().endswith("csv"):
        table_dct = pd.read_csv(uploaded_file.file).to_dict(orient="list")
    elif usr_file_name.lower().endswith("tsv"):
        table_dct = pd.read_csv(uploaded_file.file, sep="\t").to_dict(orient="list")
    else:
        if usr_file_name:
            err_lst.append("File type not supported")
        else:
            err_lst.append("No file received.")
        err_lst.append("Please upload a .csv or .xlsx file.")
        table_dct = {}
    table_dct = keep_string_only(table_dct)
    if table_dct:
        pass
    else:
        if uploaded_file.filename:
            err_lst.append(f"Can not read the uploaded file: {uploaded_file.filename}.")

    return table_dct, err_lst


def get_file_type(file_type: FileType) -> str:
    if file_type == FileType.xlsx:
        file_type = "xlsx"
    elif file_type == FileType.csv:
        file_type = "csv"
    else:
        file_type = "xlsx"

    return file_type


def get_output_name(tool: str = "", file_type: str = "xlsx") -> str:
    if tool:
        output_name = f"LipidLynxX-{tool}_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{file_type}"
    else:
        output_name = (
            f"LipidLynxX_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.{file_type}"
        )

    return output_name


def table2html(converter_data: ConverterExportData):
    data = converter_data.data
    converted_html = ""
    not_converted_html = ""
    if data:
        not_converted_dct = {}
        df_lst = []
        for k in data:
            print(type(data[k]))
            if isinstance(data[k], ConvertedListData):
                k_pairs = data[k].converted
                if k_pairs and isinstance(k, str):
                    df_lst.append(pd.DataFrame(k_pairs, columns=[k, f"{k}_converted"]))
                k_not_converted = data[k].skipped
                if k_not_converted:
                    not_converted_dct[f"{k}_skipped"] = k_not_converted

        if df_lst:
            converted_df = pd.concat(df_lst, axis=1)
            converted_df.index = converted_df.index + 1
            converted_html = converted_df.to_html()
        if not_converted_dct:
            not_converted_df = pd.DataFrame.from_dict(
                not_converted_dct, orient="index"
            ).T
            not_converted_df.index = not_converted_df.index + 1
            not_converted_html = not_converted_df.to_html()
    return converted_html, not_converted_html


def save_table(df: pd.DataFrame, file_name: str) -> (bool, str):
    is_output = False
    abs_output_path = None
    if not df.empty:
        df.to_excel(file_name)
        is_output = True
        abs_output_path = get_abs_path(file_name)

    return is_output, abs_output_path
