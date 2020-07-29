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

import logging

from lynx.utils.cfg_reader import app_log_level, cli_log_level, get_log_level


date_fmt = "%b-%d@%H:%M:%S"
log_fmt = "%(asctime)s[%(levelname)-5s] %(message)s"
app_logger = logging.getLogger("app_log")
app_logger.setLevel(app_log_level)

if not app_logger.handlers:
    app_log_handler = logging.StreamHandler()
    app_log_handler.setFormatter(fmt=logging.Formatter(log_fmt, datefmt=date_fmt))
    app_logger.info("LipidLynxX Log started ...")
    app_logger.addHandler(app_log_handler)

cli_logger = logging.getLogger("log")
cli_logger.setLevel(cli_log_level)

if not cli_logger.handlers:
    cli_log_handler = logging.StreamHandler()
    cli_log_handler.setFormatter(fmt=logging.Formatter(log_fmt, datefmt=date_fmt))
    cli_logger.info("LipidLynxX CLI Log started ...")
    cli_logger.addHandler(cli_log_handler)


def create_log(log_level="ERROR", name="log"):

    logger = logging.getLogger(name)
    logger.setLevel(get_log_level(log_level))

    if not logger.handlers:
        log_console_handler = logging.StreamHandler()
        log_console_handler.setFormatter(
            fmt=logging.Formatter(log_fmt, datefmt=date_fmt)
        )
        logger.info(f"Log started with level {log_level}...")
        logger.addHandler(log_console_handler)

    return logger
