# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os
from typing import List

from jsonschema import Draft7Validator

from ..models.log import logger


def get_abs_path(file_path: str) -> str:
    """
    check and get absolute file path from given file
    Args:
        file_path: The relative path to a file

    Returns:
        abs_path: the absolute path of input file

    """

    abs_path = ""

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
        raise FileNotFoundError(f"Can not find file: {file_path}")

    return abs_path


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


def seg_to_str(in_list: List[str], sep: str = ",") -> str:
    """
    Combine multiple segments into one str without space and separator in the end
    Args:
        in_list: input list to be joined.
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


def js_reader(file: str) -> dict:
    file = get_abs_path(file)
    if file.lower().endswith(".json"):
        with open(file) as file_obj:
            js_obj = json.load(file_obj)
            return js_obj
    else:
        raise IOError(f"Input file: {file} is not json file")


def check_json(validator: Draft7Validator, json_obj: json) -> bool:
    is_valid = False
    if validator.is_valid(json_obj):
        logger.debug(f"JSON Schema check PASSED.")
        is_valid = True
    else:
        for e in validator.iter_errors(json_obj):
            logger.error(e)
        raise Exception(f"JSON Schema check FAILED. Json string: {json_obj}")
    return is_valid
