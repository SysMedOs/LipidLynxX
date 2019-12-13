# -*- coding: utf-8 -*-
#
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
import os
import time

import pandas as pd

from ..config import app_cfg_dct
from ..controllers.general_functions import get_abs_path


def get_table(file_path: str) -> dict:

    abs_path = get_abs_path(file_path)

    if abs_path.lower().endswith(".csv"):
        df = pd.read_csv(abs_path)
    elif abs_path.lower().endswith(".tsv"):
        df = pd.read_csv(abs_path, sep="\t")
    elif abs_path.lower().endswith(".xlsx") or abs_path.lower().endswith(".xls"):
        df = pd.read_excel(abs_path)
    else:
        df = pd.DataFrame()

    if not df.empty:
        dct = df.to_dict(orient="list")
    else:
        dct = {}

    return dct


def save_table(df: pd.DataFrame, file_name: str):
    is_output = False
    abs_output_path = None
    if not df.empty:
        df.to_excel(file_name)
        is_output = True
        abs_output_path = get_abs_path(file_name)

    return is_output, abs_output_path


def clean_dct(dct: dict) -> dict:

    if dct:
        for k in dct:
            v = dct[k]
            if isinstance(v, str):
                dct[k] = [v]
            elif isinstance(v, list):
                dct[k] = list(filter(None, v))
            else:
                pass
    else:
        pass

    return dct


def create_output(data: dict) -> str:
    output_name = None
    converted_df = pd.DataFrame()
    not_converted_df = pd.DataFrame()
    if data:
        not_converted_dct = {}
        df_lst = []
        for k in data:
            k_pairs = data[k].get("PAIR", None)
            k_not_converted = data[k].get("NOT_CONVERTED", None)
            if k_pairs and isinstance(k, str):
                df_lst.append(pd.DataFrame(k_pairs, columns=[k, f"{k}_CONVERTED"]))

            if k_not_converted:
                not_converted_dct[f"{k}_NOT_CONVERTED"] = k_not_converted

        if df_lst:
            converted_df = pd.concat(df_lst, axis=1)

        if not_converted_dct:
            not_converted_df = pd.DataFrame.from_dict(
                not_converted_dct, orient="index"
            ).T

        if not converted_df.empty:
            output_name = f"LipidLynx_Output_{int(time.time())}.xlsx"
            output_path = os.path.join(app_cfg_dct["ABS_DOWNLOAD_PATH"], output_name)
            with pd.ExcelWriter(output_path, engine="openpyxl") as output_writer:
                converted_df.to_excel(output_writer, sheet_name="Converted")
                if not not_converted_df.empty:
                    not_converted_df.to_excel(output_writer, sheet_name="Not_Converted")

    else:
        pass

    return output_name
