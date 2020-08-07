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

from datetime import datetime
import json
from io import BytesIO
import os
import re
from pathlib import Path
import time
from typing import List, Union

from fastapi import File
from natsort import natsorted
import pandas as pd

from lynx.models.api_models import ConverterExportData, ConvertedListData, FileType
from lynx.utils.basics import get_abs_path
from lynx.utils.log import app_logger
from lynx.utils.toolbox import keep_string_only


def get_json(file: str) -> dict:
    file = get_abs_path(file)
    if file.lower().endswith(".json"):
        with open(file) as file_obj:
            js_obj = json.load(file_obj)
            return js_obj
    else:
        raise IOError(f"Input file: {file} is not json file")


def load_folder(folder: str, file_type: str = "", logger=app_logger) -> List[str]:
    """
     Load all files under given folder, optional with selected file suffix
     Args:
         folder: path of the folder.
         file_type: type of the file, default value is "" for no file type filter.
         logger: logger for logging.

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
    data: dict,
    output_name: Union[str, Path] = None,
    file_type: str = ".xlsx",
    converted_only: bool = False,
) -> Union[BytesIO, str]:
    file_info = None
    converted_df = pd.DataFrame()
    not_converted_df = pd.DataFrame()
    if data and not converted_only:
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
    elif data and converted_only:
        converted_df = pd.DataFrame(data)
    else:
        pass

    if not converted_df.empty:
        if output_name:
            try:
                err_msg = None
                if isinstance(output_name, Path):
                    output_name = output_name.as_posix()
                elif isinstance(output_name, str):
                    pass
                else:
                    err_msg = (
                        f"[Type error] Cannot create file: {output_name} as output."
                    )
                if output_name.lower().endswith("csv"):
                    converted_df.to_csv(output_name)
                else:
                    converted_df.to_excel(
                        output_name, sheet_name="converted", index=False
                    )
                if err_msg:
                    file_info = err_msg
                else:
                    file_info = get_abs_path(output_name)
            except IOError:
                file_info = f"[IO error] Cannot create file: {output_name} as output."
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
    data: dict, output_name: Union[str, Path] = None
) -> Union[BytesIO, str]:
    file_info = None
    is_file_name = False
    if data:
        if output_name:
            try:
                err_msg = None
                if isinstance(output_name, Path):
                    output_name = output_name.as_posix()
                elif isinstance(output_name, str):
                    pass
                else:
                    err_msg = (
                        f"[Type error] Cannot create file: {output_name} as output."
                    )
                if output_name.lower().endswith(".xlsx"):
                    table_writer = pd.ExcelWriter(output_name, engine="openpyxl")
                    is_file_name = True
                    file_info = output_name
                else:
                    err_msg = (
                        f"[Type error] Cannot create file: {output_name} as output."
                    )
                if err_msg:
                    file_info = err_msg
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
            pd.DataFrame.from_dict(skipped_dct, orient="index").T.sort_index().to_excel(
                table_writer, sheet_name="skipped", index=False
            )
        else:
            pass

        table_writer.save()
        if is_file_name:
            pass
        else:
            if isinstance(file_info, BytesIO):
                file_info.seek(0)
    else:
        file_info = None
    return file_info


def create_linker_output(
    data: dict,
    output_name: Union[str, Path] = None,
    file_type: str = ".xlsx",
    export_url: bool = True,
) -> Union[BytesIO, str]:
    file_info = None
    linked_df = pd.DataFrame()
    sum_linked_resources = {}
    if data:
        idx = 1
        for lipid_name in data:
            lipid_resources = {}
            if isinstance(data[lipid_name], dict):
                lipid_resources["Input_Name"] = data[lipid_name].get("lipid_name", "")
                lipid_resources["Shorthand_Notation"] = data[lipid_name].get(
                    "shorthand_name", ""
                )
                lipid_resources["LipidLynxX"] = data[lipid_name].get("lynx_name", "")
                lipid_resources["BioPAN"] = data[lipid_name].get("biopan_name", "")
                resource_data = data[lipid_name].get("resource_data", {})
                for db_group in resource_data:
                    db_group_resources = resource_data[db_group]
                    for db in db_group_resources:
                        db_resources = db_group_resources.get(db)
                        if db_resources and isinstance(db_resources, dict):
                            if len(list(db_resources.keys())) < 2:
                                lipid_resources[db] = ";".join(
                                    list(db_resources.keys())
                                )
                                lipid_resources[f"Link_{db}"] = ";".join(
                                    [db_resources.get(i) for i in db_resources]
                                )
                            else:
                                lipid_resources[db] = json.dumps(
                                    list(db_resources.keys())
                                )
                                lipid_resources[f"Link_{db}"] = json.dumps(
                                    [db_resources.get(i) for i in db_resources]
                                )
                        else:
                            lipid_resources[db] = ""
            sum_linked_resources[idx] = lipid_resources
            idx += 1

    default_col = ["Input_Name", "Shorthand_Notation", "LipidLynxX", "BioPAN"]
    if sum_linked_resources:
        sum_df = pd.DataFrame(data=sum_linked_resources).T
        sum_df_columns = sum_df.columns.tolist()
        link_cols = []
        for col in sum_df_columns:
            if col.startswith("Link_"):
                sum_df_columns.remove(col)
                if export_url:
                    link_cols.append(col)
            elif col in default_col:
                sum_df_columns.remove(col)
        sum_df_columns = default_col + natsorted(sum_df_columns) + natsorted(link_cols)
        linked_df = pd.DataFrame(sum_df, columns=sum_df_columns)

    if not linked_df.empty:
        if output_name:
            try:
                err_msg = None
                if isinstance(output_name, Path):
                    output_name = output_name.as_posix()
                elif isinstance(output_name, str):
                    pass
                else:
                    err_msg = (
                        f"[Type error] Cannot create file: {output_name} as output."
                    )
                if output_name.lower().endswith("csv"):
                    linked_df.to_csv(output_name)
                else:
                    linked_df.to_excel(
                        output_name, sheet_name="LinkerResult", index=False
                    )
                if err_msg:
                    file_info = err_msg
                else:
                    file_info = get_abs_path(output_name)
            except IOError:
                file_info = f"[IO error] Cannot create file: {output_name} as output."
        else:
            file_info = BytesIO()
            if file_type.lower().endswith("csv"):
                file_info.write(linked_df.to_csv().encode("utf-8"))

            else:
                output_writer = pd.ExcelWriter(
                    file_info, engine="openpyxl"
                )  # write to BytesIO instead of file path
                linked_df.to_excel(
                    output_writer, sheet_name="LinkerResult", index=False
                )
                output_writer.save()
            file_info.seek(0)

    return file_info


def get_table(uploaded_file: File, err_lst: list, logger=app_logger) -> (dict, list):
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
    table_dct = keep_string_only(table_dct, logger=logger)
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


def clean_temp_folder(
    temp_dir: str = r"lynx/temp", max_days: float = 7.0, max_files: int = 99
) -> list:
    current_time = time.time()
    earliest_unix_time = current_time - max_days * 86400  # 24 * 60 * 60 == 86400
    removed_files = []
    temp_files_lst = os.listdir(temp_dir)
    file_suffix_rgx = re.compile(r"^(.*)(\.)(csv|xlsx?)$", re.IGNORECASE)
    temp_file_path_lst = [
        os.path.join(temp_dir, f) for f in temp_files_lst if file_suffix_rgx.match(f)
    ]
    temp_file_ctime_lst = [os.path.getctime(p) for p in temp_file_path_lst]
    temp_file_info_lst = list(zip(temp_file_ctime_lst, temp_file_path_lst))
    if len(temp_file_info_lst) > max_files:
        temp_file_info_lst.sort(key=lambda tup: tup[0], reverse=True)
        removed_files = temp_file_info_lst[max_files:]
    for temp_file_info in temp_file_info_lst:
        if temp_file_info[0] < earliest_unix_time:
            removed_files.append(temp_file_info)
    for temp in removed_files:
        if os.path.isfile(temp[1]):
            os.remove(temp[1])

    return removed_files
