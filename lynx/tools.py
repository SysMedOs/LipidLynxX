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

from pathlib import Path
from typing import List

import click_spinner
import typer

from lynx.controllers.converter import Converter
from lynx.controllers.equalizer import Equalizer
from lynx.models.api_models import StyleType
from lynx.utils.cli_utils import cli_get_table, cli_save_output
from lynx.utils.file_handler import (
    create_converter_output,
    create_equalizer_output,
    get_output_name,
)
from lynx.utils.log import cli_logger
from lynx.utils.cfg_reader import lynx_version
from lynx.utils.toolbox import get_levels, get_style_level


cli_app = typer.Typer(
    help=f"LipidLynxX CLI tools version {lynx_version}. Developed by developed by team SysMedOs @ University of Leipzig"
)


@cli_app.command()
def convert_lipid(
    lipid: str = typer.Argument(None),
    style: StyleType = typer.Option(
        "LipidLynxX",
        "--style",
        "-s",
        help="The export style, choose from LipidLynxX and COMP_DB. Set to LipidLynxX by default.",
    ),
    level: str = typer.Option(
        "MAX",
        "--level",
        "-l",
        help="Level of lipid identifier information for the output. Set to MAX by default.",
    ),
):
    """
    Convert one LIPID name into supported levels and export to supported style

    LIPID: lipid abbreviation

    --style: Export style. LipidLynxX or COMP_DB. Default value: LipidLynxX

    --level: LipidLynxX lipid information levels. e.g. B0, D1, S4.1 or MAX for Maximum level. Default value: MAX

    e.g. convert PLPC --style LipidLynxX --level S1

    """
    style, level = get_style_level(style, level)
    lynx_converter = Converter(style=style)
    converted_results = lynx_converter.convert_str(input_str=lipid, level=level).get(
        "output", []
    )
    if isinstance(converted_results, list) and len(converted_results) > 0:
        converted_name = converted_results[0]
    else:
        converted_name = f"Failed to convert: {lipid}"
    typer.echo(
        typer.style(
            f"Convert {lipid} into {style} style @ {level} level:", fg=typer.colors.CYAN
        )
    )
    typer.echo(converted_name)


@cli_app.command()
def convert(
    file: Path = typer.Argument(None),
    output_file: Path = typer.Option(
        None,
        "--output",
        "-o",
        help="Path for the export file. Accept .csv / .xlsx file.",
    ),
    style: StyleType = typer.Option(
        "LipidLynxX",
        "--style",
        "-s",
        help="The export style, choose from LipidLynxX and COMP_DB. Set to LipidLynxX by default.",
    ),
    level: str = typer.Option(
        "MAX",
        "--level",
        "-l",
        help="Level of lipid identifier information for the output. Set to MAX by default.",
    ),
):
    """
    Convert one .csv / .xlsx FILE containing lipid names into supported levels
    and export to supported style as .csv / .xlsx file.
    """

    table_dct = cli_get_table(file)
    style, level = get_style_level(style, level)
    lynx_converter = Converter(style=style)
    typer.echo(
        typer.style(
            f"Convert lipid names into {style} style @ {level} level.",
            fg=typer.colors.CYAN,
        )
    )
    typer.echo(f"Processing file: {file.name} ...")
    with click_spinner.spinner():
        converted_dct = lynx_converter.convert_dict(table_dct, level=level)
    if output_file:
        pass
    else:
        input_folder = file.parent
        output_file = Path.joinpath(input_folder, get_output_name("Converter", "xlsx"))
    typer.echo(f"Generating output file...")
    with click_spinner.spinner():
        output_info = create_converter_output(converted_dct, output_name=output_file)
    cli_save_output(output_info, output_file)


@cli_app.command()
def equalize(
    file: Path = typer.Argument(None),
    output_file: Path = typer.Option(
        None,
        "--output",
        "-o",
        help="Path for the export file. Accept .xlsx file only.",
    ),
    level: List[str] = typer.Option(
        ["B1"],
        "--level",
        "-l",
        help="Level of lipid identifier information for the output. Set to B1 by default.",
    ),
):
    """
    Equalize one .csv / .xlsx FILE containing lipid names into supported levels
    and export to supported style as .csv / .xlsx file.
    """
    table_dct = cli_get_table(file)
    levels = get_levels(level)
    typer.echo(
        typer.style(f"Equalize lipid names on {levels} level.", fg=typer.colors.CYAN,)
    )
    typer.echo(f"Processing file: {file.name} ...")
    with click_spinner.spinner():
        equalizer = Equalizer(table_dct, level=levels)
        equalizer_data = equalizer.cross_match()
    if output_file:
        pass
    else:
        input_folder = file.parent
        output_file = Path.joinpath(input_folder, get_output_name("Equalizer", "xlsx"))
    typer.echo(f"Generating output file...")
    with click_spinner.spinner():
        output_info = create_equalizer_output(
            equalizer_data.dict().get("data"), output_name=output_file
        )
    cli_save_output(output_info, output_file)


if __name__ == "__main__":
    cli_app()
