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

import os

import lynx.utils


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
