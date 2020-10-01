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

import asyncio
import json
import os
import re
from pathlib import Path
from typing import List, Union

import click_spinner
import typer

from lynx.utils.log import cli_logger

from lynx.controllers.converter import Converter
from lynx.controllers.equalizer import Equalizer
from lynx.controllers.linker import get_cross_links, get_lmsd_name, get_swiss_name
from lynx.models.api_models import StyleType
from lynx.utils.cli_utils import cli_get_table, cli_save_output
from lynx.utils.file_handler import (
    create_converter_output,
    create_equalizer_output,
    get_output_name,
)

from lynx.utils.cfg_reader import app_cfg_info, lynx_version
from lynx.utils.params_loader import build_input_rules, build_output_rules
from lynx.utils.toolbox import get_levels, get_style_level


cli_app = typer.Typer(
    help=f"LipidLynxX CLI tools version {lynx_version}. "
    f"Developed by team SysMedOs @ University of Leipzig "
    f"Please cite our publication in an appropriate form. "
    f"LipidLynxX preprint on bioRxiv.org. Zhixu Ni, Maria Fedorova. "
    f"'LipidLynxX: a data transfer hub to support integration of large scale lipidomics datasets' "
    f"DOI: 10.1101/2020.04.09.033894"
)

default_input_rules = build_input_rules(app_cfg_info["input_rules"], cli_logger)
default_output_rules = build_output_rules(app_cfg_info["output_rules"], cli_logger)


@cli_app.command(name="convert-lipid")
def convert_lipid(
    lipid: str = typer.Argument(None),
    style: StyleType = typer.Option(
        "LipidLynxX",
        "--style",
        "-s",
        help="The export style, choose from LipidLynxX, COMP_DB, and ShorthandNotation. Set to LipidLynxX by default.",
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

    --style: Export style. LipidLynxX, COMP_DB, or ShorthandNotation. Default value: LipidLynxX

    --level: LipidLynxX lipid information levels. e.g. B0, D1, S4.1 or MAX for Maximum level. Default value: MAX

    e.g. convert "PLPC" --style LipidLynxX --level S1

    """
    style, level = get_style_level(style, level)
    lynx_converter = Converter(
        style=style,
        input_rules=default_input_rules,
        output_rules=default_output_rules,
        logger=cli_logger,
    )
    if lipid:
        converted_name = lynx_converter.convert_str(input_str=lipid, level=level).output
        if isinstance(converted_name, str) and len(converted_name) > 0:
            pass
        else:
            converted_name = f"Failed to convert: {lipid}"
        typer.secho(
            f"Convert {lipid} into {style} style @ {level} level:",
            fg=typer.colors.CYAN,
        )
        typer.echo(converted_name)
    else:
        converted_name = ""
        typer.secho(
            'Please input a lipid name. e.g. "PLPC".', fg=typer.colors.YELLOW,
        )
        typer.echo(convert_lipid.__doc__)

    return converted_name


@cli_app.command(name="convert-lipids")
def convert_lipids(
    lipids: str = typer.Argument(None),
    style: StyleType = typer.Option(
        "LipidLynxX",
        "--style",
        "-s",
        help="The export style, choose from LipidLynxX, COMP_DB, and ShorthandNotation. Set to LipidLynxX by default.",
    ),
    level: str = typer.Option(
        "MAX",
        "--level",
        "-l",
        help="Level of lipid identifier information for the output. Set to MAX by default.",
    ),
):
    """
    Convert LIPIDS into supported levels and export to supported style

    LIPIDS: lipid abbreviation

    --style: Export style. LipidLynxX, COMP_DB, or ShorthandNotation. Default value: LipidLynxX

    --level: LipidLynxX lipid information levels. e.g. B0, D1, S4.1 or MAX for Maximum level. Default value: MAX

    e.g.
    python cli_lynx.py convert_lipids '["PLPC","DG(O-16:0/18:1)"]' --style LipidLynxX --level B1
    python cli_lynx.py convert_lipids '{"Source_1":["DHA","PLPC","PI 16:0-20:4","UNKNOWN_LIPID_1"],"Source_2":["FA20:4","PC 16:0_18:2","PI(16:0/20:4)","UNKNOWN_LIPID_2"]}' --style LipidLynxX

    """
    if lipids:
        if '"' in lipids:
            js_obj = None
            try:
                js_obj = json.loads(lipids)
            except (json.JSONDecodeError, json.decoder.JSONDecodeError):
                typer.echo(
                    typer.style(
                        f"Currently support JSON list/array dict/object object only",
                        fg=typer.colors.YELLOW,
                    )
                )
            if isinstance(js_obj, str):
                convert_lipid(lipid=js_obj, style=style, level=level)
            else:
                style, level = get_style_level(style, level)
                lynx_converter = Converter(
                    style=style,
                    input_rules=default_input_rules,
                    output_rules=default_output_rules,
                    logger=cli_logger,
                )
                if isinstance(js_obj, list):
                    converted_obj = lynx_converter.convert_list(
                        input_list=js_obj, level=level
                    )
                    converted_names = converted_obj.output
                    typer.echo(
                        typer.style(
                            f"Convert {js_obj} into {style} style @ {level} level:",
                            fg=typer.colors.CYAN,
                        )
                    )
                    typer.echo(json.dumps(converted_names))
                    skipped_names = converted_obj.skipped
                    if skipped_names:
                        typer.secho(
                            f"Skipped for conversion: {skipped_names}",
                            fg=typer.colors.YELLOW,
                        )
                elif isinstance(js_obj, dict):
                    converted_dct = lynx_converter.convert_dict(
                        input_dct=js_obj, level=level
                    )
                    typer.echo(
                        typer.style(
                            f"Convert {js_obj} into {style} style @ {level} level:",
                            fg=typer.colors.CYAN,
                        )
                    )
                    export_dct = {}
                    for k in converted_dct:
                        export_dct[k] = converted_dct[k].dict()
                    typer.echo(json.dumps(export_dct))
                else:
                    typer.echo(
                        typer.style(
                            f"Currently support JSON list/array dict/object object only",
                            fg=typer.colors.YELLOW,
                        )
                    )
        else:
            convert_lipid(lipid=lipids, style=style, level=level)
    else:
        typer.secho(
            'Please input a lipid name. e.g. "PLPC" or a file path. Type --help to see instructions.',
            fg=typer.colors.YELLOW,
        )
        typer.echo(convert_lipids.__doc__)


@cli_app.command(name="convert-file")
def convert_file(
    file: Path = typer.Argument(None),
    column: str = typer.Option(
        None, "--column", "-c", help="name of the column that contains lipid notations",
    ),
    output_file: Path = typer.Option(
        None,
        "--output",
        "-o",
        help="Path for the export file. Accept .txt/ .csv / .xlsx file.",
    ),
    style: StyleType = typer.Option(
        "LipidLynxX",
        "--style",
        "-s",
        help="The export style, choose from LipidLynxX, COMP_DB, and ShorthandNotation. Set to LipidLynxX by default.",
    ),
    level: str = typer.Option(
        "MAX",
        "--level",
        "-l",
        help="Level of lipid identifier information for the output. Set to MAX by default.",
    ),
):
    """
    Convert one .txt/ .csv / .xlsx FILE containing lipid names into supported levels
    and export to supported style as .csv / .xlsx file.
    """
    if column:
        lipid_col_name = column
        use_one_col = True
        converted_only = True
    else:
        lipid_col_name = None
        use_one_col = False
        converted_only = False

    raw_table_dct = cli_get_table(file)
    if lipid_col_name:
        table_dct = {lipid_col_name: raw_table_dct.get(lipid_col_name)}
        if table_dct.get(lipid_col_name):
            pass
        else:
            typer.echo(
                typer.style(
                    f"Can NOT find column '{lipid_col_name}' in file: {file}",
                    fg=typer.colors.YELLOW,
                )
            )
    else:
        table_dct = raw_table_dct
    style, level = get_style_level(style, level)
    lynx_converter = Converter(
        style=style,
        input_rules=default_input_rules,
        output_rules=default_output_rules,
        logger=cli_logger,
    )
    typer.echo(
        typer.style(
            f"Convert lipid names into {style} style @ {level} level.",
            fg=typer.colors.CYAN,
        )
    )
    typer.echo(f"Processing file: {file.name} ...")

    with click_spinner.spinner():
        if use_one_col:
            lipid_list = table_dct.get(lipid_col_name)
            converted_obj = lynx_converter.convert_list(
                input_list=lipid_list, level=level, default_na="UNPROCESSED"
            )
            converted_names = converted_obj.output
            converted_dct = {f"Converted_{lipid_col_name}": converted_names}
            converted_dct.update(raw_table_dct)
        else:
            converted_dct = lynx_converter.convert_dict(table_dct, level=level)
    if output_file:
        if isinstance(output_file, Path):
            pass
        else:
            output_file = Path(output_file)
    else:
        input_folder = file.parent
        output_file = Path.joinpath(input_folder, get_output_name("Converter", "xlsx"))
    typer.echo(f"Generating output file...")
    with click_spinner.spinner():
        output_info = create_converter_output(
            converted_dct, output_name=output_file, converted_only=converted_only
        )

    cli_save_output(output_info, output_file)


@cli_app.command(name="convert")
def convert(
    source: str = typer.Argument(None),
    column: str = typer.Option(
        None, "--column", "-c", help="name of the column that contains lipid notations",
    ),
    style: StyleType = typer.Option(
        "LipidLynxX",
        "--style",
        "-s",
        help="The export style, choose from LipidLynxX, COMP_DB, and ShorthandNotation. Set to LipidLynxX by default.",
    ),
    level: str = typer.Option(
        "MAX",
        "--level",
        "-l",
        help="Level of lipid identifier information for the output. Set to MAX by default.",
    ),
    output_file: Path = typer.Option(
        None,
        "--output",
        "-o",
        help=(
            "Path for the export file. Accept .csv / .xlsx file. "
            "This option is valid only when an input file is provided as SOURCE argument"
        ),
    ),
):
    """
    Convert SOURCE lipid names into supported levels.

    Support SOURCE input as JSON string or file path of .txt/ .csv/ .xlsx file.

    Use option --output to export a .csv / .xlsx file from input file.

    e.g.

    *convert one name:

    python cli_lynx.py convert "DG(O-16:0/18:1)" -s COMP_DB

    * convert a List of names:

    python cli_lynx.py convert '["PLPC","DG(O-16:0/18:1)"]' --style LipidLynxX --level B1

    * convert a Dict of names in List:

    python cli_lynx.py convert '{"Source_1":["DHA","PLPC","PI 16:0-20:4","UNKNOWN_LIPID_1"],"Source_2":["FA20:4","PC 16:0_18:2","PI(16:0/20:4)","UNKNOWN_LIPID_2"]}' --style LipidLynxX

    * convert an input file and export to an output file:

    python cli_lynx.py convert path/to/input_file.csv --output path/to/output_file.xlsx -s COMP_DB

    """
    if source:

        if isinstance(source, Path):
            convert_file(file=source, output_file=output_file, style=style, level=level)
        elif isinstance(source, str):
            if re.match(r".+(\.csv|\.txt|\.xlsx)$", source, re.IGNORECASE):
                if os.path.isfile(source):
                    convert_file(
                        file=Path(source),
                        column=column,
                        output_file=output_file,
                        style=style,
                        level=level,
                    )
            else:
                convert_lipids(lipids=source, style=style, level=level)
        else:
            typer.secho(
                'Please input a lipid name. e.g. "PLPC" or a file path. Type --help to see instructions.',
                fg=typer.colors.YELLOW,
            )
            typer.echo(convert.__doc__)
    else:
        typer.secho(
            'Please input a lipid name. e.g. "PLPC" or a file path. Type --help to see instructions.',
            fg=typer.colors.YELLOW,
        )
        typer.echo(convert.__doc__)


@cli_app.command(name="equalize")
@cli_app.command(name="equalize-file")
def equalize(
    file: Path = typer.Argument(None),
    output_file: Path = typer.Option(
        None,
        "--output",
        "-o",
        help="Path for the export file. Accept .xlsx file only.",
    ),
    level: List[str] = typer.Option(
        ["B2"],
        "--level",
        "-l",
        help="Level of lipid identifier information for the output. Set to B2 by default.",
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
        equalizer = Equalizer(
            table_dct,
            level=levels,
            input_rules=default_input_rules,
            output_rules=default_output_rules,
            logger=cli_logger,
        )
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


@cli_app.command(name="link-lipid")
def link_lipid(
    lipid: str = typer.Argument(None),
    formatted: bool = typer.Option(False, "--formatted", "-f"),
    url: bool = typer.Option(False, "--url", "-u"),
):
    """
    Link one LIPID name into public available databases and resources

    LIPID: lipid abbreviation

    --formatted: Group results by types of resources

    --url: export URL of the resources

    e.g. link-lipid "PLPC" --formatted --url

    """
    if lipid:
        if re.match(r"^LM\w\w\d{8}$", lipid, re.IGNORECASE):
            safe_lipid_name = asyncio.run(get_lmsd_name(lipid))
        elif re.match(r"^SLM:\d{9}$", lipid, re.IGNORECASE):
            safe_lipid_name = asyncio.run(get_swiss_name(lipid))
        else:
            safe_lipid_name = lipid
        search_name = convert_lipid(
            safe_lipid_name, style=StyleType("BracketsShorthand"), level="MAX"
        )
        result = asyncio.run(
            get_cross_links(lipid_name=search_name, export_url=url, formatted=formatted)
        )
        typer.secho(json.dumps(result))
    else:
        typer.secho(
            'Please input a lipid name. e.g. "PLPC". Type --help to see instructions.',
            fg=typer.colors.YELLOW,
        )


if __name__ == "__main__":
    cli_app()
