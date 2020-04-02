# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import logging

log_level = logging.DEBUG
date_fmt = "%b-%d@%H:%M:%S"
log_fmt = "%(asctime)s[%(levelname)-5s] %(message)s"
logger = logging.getLogger("log")
logger.setLevel(log_level)

if not logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(fmt=logging.Formatter(log_fmt, datefmt=date_fmt))
    logger.info("Log started ...")
    logger.addHandler(console_handler)
