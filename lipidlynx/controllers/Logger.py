# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
import logging

logger = logging.getLogger("log")
log_level = logging.DEBUG
logger.setLevel(log_level)
if not logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.WARNING)
    formatter = logging.basicConfig(
        format="%(asctime)s[%(levelname)-5s] %(message)s",
        datefmt="%b-%d@%H:%M:%S",
        level=log_level,
    )
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
