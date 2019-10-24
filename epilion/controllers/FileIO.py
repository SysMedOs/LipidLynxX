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
import pandas as pd

from epilion.controllers.GeneralFunctions import get_abs_path


def get_table(file_path):

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
