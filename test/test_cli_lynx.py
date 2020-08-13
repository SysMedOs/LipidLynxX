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
import sys

from typer.testing import CliRunner

lynx_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lynx_path + "/../")

from cli_lynx import cli_app
from lynx.utils.file_handler import get_abs_path


runner = CliRunner()


def test_convert_lipid():
    result = runner.invoke(cli_app, ["convert-lipid", "PC 16:0/18:2"])
    # assert result.exit_code == 0
    print(result.stdout)
    assert "PC(16:0/18:2)" in result.stdout


def test_convert():
    try:
        test_output = get_abs_path(r"test/test_output/test_convert.xlsx")
    except FileNotFoundError:
        if os.path.isdir(r"test/test_output"):
            test_output = os.path.join(r"test/test_output", "test_convert.xlsx")
            test_output = os.path.abspath(test_output)
        else:
            test_output = r"test/test_output/test_convert.xlsx"
    if os.path.isfile(test_output):
        os.remove(test_output)
    result = runner.invoke(
        cli_app,
        [
            "convert",
            r"doc/sample_data/input/LipidLynxX_test.csv",
            "--output",
            test_output,
        ],
    )
    print(result.stdout)
    # assert result.exit_code == 0
    # assert "Save output as: test/test_output/test_convert_cli.xlsx" in result.stdout
    assert os.path.isfile(test_output)


def test_equalize():
    try:
        test_output = get_abs_path(r"test/test_output/test_equalize.xlsx")
    except FileNotFoundError:
        if os.path.isdir(r"test/test_output"):
            test_output = os.path.join(r"test/test_output", "test_equalize.xlsx")
            test_output = os.path.abspath(test_output)
        else:
            test_output = r"test/test_output/test_equalize.xlsx"
    if os.path.isfile(test_output):
        os.remove(test_output)
    result = runner.invoke(
        cli_app,
        [
            "equalize",
            r"doc/sample_data/input/LipidLynxX_test.csv",
            "--output",
            test_output,
        ],
    )
    print(result.stdout)
    # assert result.exit_code == 0
    # assert "Save output as: test/test_output/test_equalize_cli.xlsx" in result.stdout
    assert os.path.isfile(test_output)
