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
        test_output = os.path.abspath(
            get_abs_path(r"test/test_output/test_convert_cli.xlsx")
        )
    except FileNotFoundError:
        if os.path.isdir(r"test/test_output"):
            test_output = os.path.join(r"test/test_output", "test_convert_cli.xlsx")
            test_output = os.path.abspath(test_output)
        else:
            test_output = r"test/test_output/test_convert_cli.xlsx"
    test_input_file = get_abs_path(r"doc/sample_data/input/LipidLynxX_test.csv")
    if os.path.isfile(test_input_file):
        import pandas as pd

        df = pd.read_csv(test_input_file)
        print(df.head())
    else:
        print("can not load input file")
    if os.path.isfile(test_output):
        os.remove(test_output)
    result = runner.invoke(
        cli_app, ["convert", test_input_file, "--output", test_output],
    )
    cli_output_lst = result.stdout.strip("\n").split("\n")
    print(cli_output_lst)
    print(os.path.abspath(test_output))
    if os.path.isfile(test_output):
        print("output created")
        import pandas as pd

        df = pd.read_excel(test_output)
        print(df.head())
    assert os.path.isfile(test_output) is True
    assert f"Save output as: {os.path.abspath(test_output)}" in cli_output_lst


def test_equalize():
    try:
        test_output = os.path.abspath(
            get_abs_path(r"test/test_output/test_equalize.xlsx")
        )
    except FileNotFoundError:
        if os.path.isdir(r"test/test_output"):
            test_output = os.path.join(r"test/test_output", "test_equalize.xlsx")
            test_output = os.path.abspath(test_output)
        else:
            test_output = os.path.abspath(r"test/test_output/test_equalize.xlsx")
    test_input_file = get_abs_path(r"doc/sample_data/input/LipidLynxX_test.csv")
    if os.path.isfile(test_output):
        os.remove(test_output)
    result = runner.invoke(
        cli_app, ["equalize", test_input_file, "--output", test_output],
    )
    cli_output_lst = result.stdout.strip("\n").split("\n")
    print(cli_output_lst)
    print(os.path.abspath(test_output))
    assert os.path.isfile(test_output) is True
    assert f"Save output as: {os.path.abspath(test_output)}" in cli_output_lst
