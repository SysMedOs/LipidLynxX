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

from lynx.models.defaults import (
    default_temp_folder,
    default_temp_max_days,
    default_temp_max_files,
)
from lynx.utils.log import app_logger
from lynx.utils.file_handler import clean_temp_folder


removed_files = clean_temp_folder(
    default_temp_folder, default_temp_max_days, default_temp_max_files
)
if removed_files:
    app_logger.info(
        f"Remove temporary output files older than {default_temp_max_days} days..."
    )
    for removed_file in removed_files:
        app_logger.info(f"File removed: {removed_file}")
