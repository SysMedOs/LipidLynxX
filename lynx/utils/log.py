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

import logging

from lynx.utils.cfg_reader import app_log_level, cli_log_level


date_fmt = "%b-%d@%H:%M:%S"
log_fmt = "%(asctime)s[%(levelname)-5s] %(message)s"
app_logger = logging.getLogger("app_log")
app_logger.setLevel(app_log_level)

if not app_logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(fmt=logging.Formatter(log_fmt, datefmt=date_fmt))
    app_logger.info("Log started ...")
    app_logger.addHandler(console_handler)

cli_logger = logging.getLogger("log")
cli_logger.setLevel(cli_log_level)

if not cli_logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(fmt=logging.Formatter(log_fmt, datefmt=date_fmt))
    cli_logger.info("Log started ...")
    cli_logger.addHandler(console_handler)
