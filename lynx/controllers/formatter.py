# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.n

import os

from natsort import natsorted
import regex as re

from lynx.models.cv import CV
from lynx.models.defaults import default_cv_file
from lynx.utils.log import logger


class Formatter(object):
    def __init__(self, cv_file: str = default_cv_file):
        if os.path.isfile(cv_file):
            pass
        else:
            cv_file = default_cv_file
        self.cv = CV(cv_file).info

    def format_mod_type(self, info: dict) -> dict:
        mod_type_lst = info.get("MOD_TYPE", [])
        mod_site_lst = info.get("MOD_SITE", [])
        mod_site_info_lst = info.get("MOD_SITE_INFO", [])
        mod_site_lst = [s.strip(" ") for s in mod_site_lst]
        mod_site_lst = [s.strip(",") for s in mod_site_lst]
        mod_site_info_lst = [si.strip(" ") for si in mod_site_info_lst]
        formatted_mod_type_lst = []
        mod_lv_dct = {}
        for mod_type in mod_type_lst:
            for alia in self.cv:
                alia_rgx = self.cv[alia].get("MATCH", None)
                matched_cv = self.cv[alia].get("CV", None)
                if alia_rgx and matched_cv is not None:
                    alia_match = alia_rgx.match(mod_type)
                    if alia_match:
                        logger.info(f"mod_type: {mod_type} identified as {matched_cv}")
                        formatted_mod_type_lst.append(matched_cv)
                        mod_lv_dct[matched_cv] = max(
                            mod_lv_dct.get(matched_cv, 0), self.cv[alia].get("LEVEL", 0)
                        )

        if (
            mod_type_lst
            and formatted_mod_type_lst
            and len(mod_type_lst) == formatted_mod_type_lst
        ):
            info["MOD_TYPE"] = formatted_mod_type_lst

        if len(formatted_mod_type_lst) == len(mod_site_lst) and len(
            formatted_mod_type_lst
        ) == len(mod_site_info_lst):
            formatted_mod_lst = zip(
                formatted_mod_type_lst, mod_site_lst, mod_site_info_lst
            )
        else:
            formatted_mod_lst = []

        mod_info_dct = {}
        if formatted_mod_lst:
            for mod_tp in formatted_mod_lst:
                mod_type = mod_tp[0]
                existed_mod_count = mod_info_dct.get(mod_type, {}).get("MOD_COUNT", 0)
                existed_mod_site_lst = mod_info_dct.get(mod_type, {}).get(
                    "MOD_SITE", []
                )
                existed_mod_site_info_lst = mod_info_dct.get(mod_type, {}).get(
                    "MOD_SITE_INFO", []
                )
                existed_mod_site_lst.append(mod_tp[1]),
                existed_mod_site_info_lst.append(mod_tp[2])
                mod_level = mod_lv_dct.get(mod_type, 0)
                mod_count = existed_mod_count + 1
                if mod_type == "DB":
                    true_site_lst = [s for s in existed_mod_site_lst if s != ""]
                    true_site_info_lst = [
                        s for s in existed_mod_site_info_lst if s != ""
                    ]
                    if (
                        len(true_site_lst) == mod_count
                        and len(true_site_info_lst) != mod_count
                    ):
                        if (
                            "Z" not in true_site_info_lst
                            and "E" not in true_site_info_lst
                        ):
                            mod_level += 0.1
                    elif (
                        len(true_site_lst) == mod_count
                        and len(true_site_info_lst) == mod_count
                    ):
                        if "Z" in true_site_info_lst or "E" in true_site_info_lst:
                            mod_level += 0.2
                    else:
                        pass
                elif mod_type != "DB" and mod_level >= 3:
                    true_site_lst = [s for s in existed_mod_site_lst if s != ""]
                    true_site_info_lst = [
                        s for s in existed_mod_site_info_lst if s != ""
                    ]
                    if (
                        len(true_site_lst) == mod_count
                        and len(true_site_info_lst) != mod_count
                    ):
                        if (
                            "S" not in true_site_info_lst
                            and "R" not in true_site_info_lst
                        ):
                            mod_level += 1
                    elif (
                        len(true_site_lst) == mod_count
                        and len(true_site_info_lst) == mod_count
                    ):
                        if "S" in true_site_info_lst or "R" in true_site_info_lst:
                            mod_level += 2
                    else:
                        pass
                updated_mod_info = {
                    "MOD_LEVEL": mod_level,
                    "MOD_COUNT": mod_count,
                    "MOD_SITE": existed_mod_site_lst,
                    "MOD_SITE_INFO": existed_mod_site_info_lst,
                }
                mod_info_dct[mod_type] = updated_mod_info

        return mod_info_dct

    def format_residue(self, info: dict) -> dict:
        residue_info_dct = {}
        link_lst = info.get("LINK", [""])
        if link_lst:
            link = link_lst[0]
        else:
            link = ""
        num_o_lst = info.get("NUM_O", ["0"])
        if num_o_lst:
            num_o = int(num_o_lst[0])
        else:
            num_o = 0

        residue_info_dct["LINK"] = link
        residue_info_dct["MOD"] = self.format_mod_type(info)
        residue_info_dct["NUM_C"] = int(info.get("NUM_C", ["0"])[0])
        residue_info_dct["NUM_DB"] = int(info.get("DB", ["0"])[0])
        residue_info_dct["NUM_O"] = num_o

        return residue_info_dct

    def format(self, info: dict) -> dict:
        formatted_info = {}

        formatted_info["RESIDUE"] = self.format_residue(info)

        return formatted_info


if __name__ == "__main__":

    pass
