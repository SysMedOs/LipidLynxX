# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os


def get_abs_path(file_path: str) -> str:

    abs_path = ""

    in_file_lst = [
        r"{path}".format(path=file_path),
        r"./{path}".format(path=file_path),
        r"../{path}".format(path=file_path),
        r"../../{path}".format(path=file_path),
        r"../../../{path}".format(path=file_path),
    ]
    for f in in_file_lst:
        if os.path.isfile(f):
            abs_path = os.path.abspath(f)
            break

    if not abs_path:
        raise FileNotFoundError(f"Can not find file: {file_path}")

    return abs_path
