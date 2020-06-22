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
import os, sys
from typer.testing import CliRunner

lynx_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lynx_path + "/../")

from cli_lynx import cli_app


runner = CliRunner()


def test_convert_lipid():
    result = runner.invoke(cli_app, ["convert-lipid", "PLPC"])
    assert result.exit_code == 0
    assert "PC(16:0/18:2)" in result.stdout


def test_convert():
    result = runner.invoke(
        cli_app,
        [
            "convert",
            r"doc/sample_data/input/LipidLynxX_test.csv",
            "--output",
            r"test/test_output/test_convert.xlsx",
        ],
    )
    assert result.exit_code == 0
    assert "Save output as: test/test_output/test_convert.xlsx" in result.stdout


def test_equalize():
    result = runner.invoke(
        cli_app,
        [
            "equalize",
            r"doc/sample_data/input/LipidLynxX_test.csv",
            "--output",
            r"test/test_output/test_equalize.xlsx",
        ],
    )
    assert result.exit_code == 0
    assert "Save output as: test/test_output/test_equalize.xlsx" in result.stdout
