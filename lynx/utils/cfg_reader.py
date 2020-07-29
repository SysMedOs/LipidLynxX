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

import configparser
import logging
import re
from typing import Dict, Union

from lynx.utils.basics import get_abs_path


def load_cfg_info(cfg_path: str = None) -> Dict[str, str]:
    cfg_dct = {}
    default_fields = [
        "api_version",
        "app_log_level",
        "app_url",
        "app_port",
        "cli_log_level",
        "controlled_vocabularies",
        "defined_alias",
        "input_rules",
        "output_rules",
        "resource_kegg",
        "resource_lion",
    ]
    config = configparser.ConfigParser()
    if cfg_path and isinstance(cfg_path, str):
        config_path = get_abs_path(cfg_path)
    else:
        try:
            config_path = get_abs_path("config.ini")
        except FileNotFoundError:
            config_path = get_abs_path("configure.ini")

    config.read(config_path)
    if config.has_section("settings"):
        user_cfg = "settings"
    elif config.has_section("default"):
        user_cfg = "default"
    else:
        user_cfg = ""
        raise ValueError(f"Cannot __load__ settings from file {config_path}")

    if len(user_cfg) > 2:
        options = config.options(user_cfg)
        for field in default_fields:
            if field in options and field in [
                "controlled_vocabularies",
                "defined_alias",
                "input_rules",
                "output_rules",
            ]:
                cfg_dct[field] = get_abs_path(config.get(user_cfg, field))
            else:
                cfg_dct[field] = config.get(user_cfg, field)

    if "app_url" not in cfg_dct:
        cfg_dct["app_url"] = "127.0.0.1"
    if "app_port" not in cfg_dct:
        cfg_dct["app_port"] = "1399"

    return cfg_dct


def get_log_level(log_level: str = "DEBUG",):
    if re.search(r"DEBUG", log_level, re.IGNORECASE):
        level = logging.DEBUG
    elif re.search(r"CRITI", log_level, re.IGNORECASE):
        level = logging.CRITICAL
    elif re.search(r"ERR", log_level, re.IGNORECASE):
        level = logging.ERROR
    elif re.search(r"INFO", log_level, re.IGNORECASE):
        level = logging.INFO
    elif re.search(r"WARN", log_level, re.IGNORECASE):
        level = logging.WARNING
    else:
        level = logging.DEBUG
    return level


def get_cli_log_settings(log_cfg: Union[str, bool] = "OFF"):
    if isinstance(log_cfg, bool):
        pass
    elif isinstance(log_cfg, str):
        if re.search(r"OFF", log_cfg, re.IGNORECASE):
            log_cfg = False
        elif re.search(r"False", log_cfg, re.IGNORECASE):
            log_cfg = False
        elif re.search(r"ON", log_cfg, re.IGNORECASE):
            log_cfg = True
        elif re.search(r"True", log_cfg, re.IGNORECASE):
            log_cfg = True
        else:
            log_cfg = False
    else:
        log_cfg = False

    return log_cfg


default_cfg_path = "/lynx/config.ini"
app_cfg_info = load_cfg_info(cfg_path=default_cfg_path)

lynx_version = "0.8.0"
api_version = app_cfg_info.get("api_version", "1.0")
app_log_level = get_log_level(app_cfg_info.get("app_log_level", "DEBUG"))
cli_log_level = get_log_level(app_cfg_info.get("cli_log_level", "ERROR"))
