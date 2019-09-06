# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os
from typing import List


def get_abs_path(file_path: str) -> str:
    """
    check and get absolute file path from given file
    Args:
        file_path: The relative path to a file

    Returns:
        abs_path: the absolute path of input file

    """

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


def seg_to_str(in_list: List[str], sep: str = ",") -> str:
    """
    Combine multiple segments into one str without space and separator in the end
    Args:
        in_list: input list to be joined
        sep: separator used to join str segments, default value is ","

    Returns:
        out_str: the joined str

    """
    in_list = filter(None, in_list)
    out_str = sep.join(in_list)
    out_str = out_str.strip(sep)
    if out_str is None:
        out_str = ""
    return out_str
