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
import os
import re
import time

from lynx.models.defaults import (
    default_temp_folder,
    default_temp_max_days,
    default_temp_max_files,
)
from lynx.utils.log import app_logger


def clean_temp_folder(
    temp_dir: str = r"lynx/temp", max_days: float = 7.0, max_files: int = 99
) -> list:
    current_time = time.time()
    earliest_unix_time = current_time - max_days * 86400  # 24 * 60 * 60 == 86400
    removed_files = []
    temp_files_lst = os.listdir(temp_dir)
    file_suffix_rgx = re.compile(r"^(.*)(\.)(csv|xlsx|txt|json?)$", re.IGNORECASE)
    temp_file_path_lst = [
        os.path.join(temp_dir, f) for f in temp_files_lst if file_suffix_rgx.match(f)
    ]
    temp_file_ctime_lst = [os.path.getctime(p) for p in temp_file_path_lst]
    temp_file_info_lst = list(zip(temp_file_ctime_lst, temp_file_path_lst))
    if len(temp_file_info_lst) > max_files:
        temp_file_info_lst.sort(key=lambda tup: tup[0], reverse=True)
        removed_files = temp_file_info_lst[max_files:]
    for temp_file_info in temp_file_info_lst:
        if temp_file_info[0] < earliest_unix_time:
            removed_files.append(temp_file_info)
    for temp in removed_files:
        if os.path.isfile(temp[1]):
            os.remove(temp[1])

    return removed_files


def remove_temp_file():
    removed_files = clean_temp_folder(
        default_temp_folder, default_temp_max_days, default_temp_max_files
    )
    if removed_files:
        app_logger.info(
            f"Remove temporary output files older than {default_temp_max_days} days..."
        )
        for removed_file in removed_files:
            app_logger.info(f"File removed: {removed_file}")
    else:
        app_logger.debug(
            f"No outdated temporary files to be removed."
        )

