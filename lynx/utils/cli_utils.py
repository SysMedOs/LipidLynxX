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

from io import BytesIO
import os
from pathlib import Path
from typing import Union

import pandas as pd
import typer


def cli_get_table(file: Union[Path, str]):
    if isinstance(file, Path) and file.is_file():
        in_file_name_low_str = file.name.lower()
    elif isinstance(file, str) and os.path.isfile(file):
        in_file_name_low_str = file.lower()
    else:
        typer.echo(f"[IO Error] Can not find file: {file}")
        raise typer.Exit(code=1)
    if in_file_name_low_str:
        if in_file_name_low_str.endswith("xlsx"):
            table_dct = pd.read_excel(file).to_dict(orient="list")
        elif in_file_name_low_str.endswith("csv"):
            table_dct = pd.read_csv(file).to_dict(orient="list")
        elif in_file_name_low_str.endswith("tsv"):
            table_dct = pd.read_csv(file, sep="\t").to_dict(orient="list")
        else:
            typer.echo("File type not supported")
            raise typer.Exit(code=1)
    else:
        typer.echo(f"[IO Error] Can not find file.")
        raise typer.Exit(code=1)

    return table_dct


def cli_save_output(output_info: Union[str, BytesIO], output_file: Path):
    if isinstance(output_info, str):
        typer.echo(
            typer.style(
                f"Save output as: {output_file.as_posix()}", fg=typer.colors.CYAN
            )
        )
        raise typer.Exit(code=0)
    else:
        typer.echo(f"Failed to generate output: {output_file.as_posix()}")
        raise typer.Exit(code=1)
