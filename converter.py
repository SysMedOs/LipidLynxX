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

from typing import Optional

import typer

from lynx.controllers.converter import Converter
from lynx.models.api_models import FileType, StyleType
from lynx.utils.toolbox import get_style_level

converter_app = typer.Typer()


@converter_app.command()
def convert(lipid: str, style: Optional[StyleType] = "LipidLynxX", level: Optional[str] = "MAX"):
    """
    Convert one LIPID name into supported levels and export to supported style

    LIPID: lipid abbreviation

    --style: Export style. LipidLynxX or COMP_DB. Default value: LipidLynxX

    --level: LipidLynxX lipid information levels. e.g. B0, D1, S4.1 or MAX for Maximum level. Default value: MAX

    e.g. convert PLPC --style LipidLynxX --level S1

    """
    style, level = get_style_level(style, level)
    lynx_converter = Converter(style=style)
    converted_results = lynx_converter.convert_str(
        input_str=lipid, level=level
    ).get("output", [])
    if isinstance(converted_results, list) and len(converted_results) > 0:
        converted_name = converted_results[0]
    else:
        converted_name = f"Failed to convert: {lipid}"
    typer.echo(f'Convert {lipid} into {style} style @ {level} level:')
    typer.echo(converted_name)


@converter_app.command()
def convert_file(input_file: str, output_file: str, file_format: FileType = FileType.xlsx):
    """
    Convert one .csv / .xlsx file of lipid names into supported levels and export to supported style
    """
    typer.echo(input_file)
    typer.echo(output_file)
    typer.echo(file_format)


if __name__ == "__main__":

    from lynx.utils.log import logger

    logger.disabled = True
    converter_app()
