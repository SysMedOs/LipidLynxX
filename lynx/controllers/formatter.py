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
from lynx.models.defaults import default_cv_file, elem_nominal_info
from lynx.utils.log import logger


class Formatter(object):
    def __init__(self, cv_file: str = default_cv_file):
        if os.path.isfile(cv_file):
            pass
        else:
            cv_file = default_cv_file
        self.alias2cv = CV(cv_file).info
        self.raw_cv = CV(cv_file).raw_cv

    @staticmethod
    def to_mass_shift(elements: dict) -> int:

        delta = 0
        for elem in elements:
            if elem in elem_nominal_info:
                delta += elem_nominal_info[elem] * elements[elem]
        return delta

    def format_mods(self, info: dict) -> dict:
        mod_type_lst = info.get("MOD_TYPE", [])
        mod_site_lst = info.get("MOD_SITE", [])
        mod_site_info_lst = info.get("MOD_SITE_INFO", [])
        mod_site_lst = [s.strip(" ") for s in mod_site_lst]
        mod_site_lst = [s.strip(",") for s in mod_site_lst]
        mod_site_info_lst = [si.strip(" ") for si in mod_site_info_lst]
        formatted_mod_type_lst = []
        mod_lv_dct = {}
        mass_shift_dct = {}
        for mod_type in mod_type_lst:
            for alia in self.alias2cv:
                alia_rgx = self.alias2cv[alia].get("MATCH", None)
                matched_cv = self.alias2cv[alia].get("CV", None)
                if alia_rgx and matched_cv is not None:
                    alia_match = alia_rgx.match(mod_type)
                    if alia_match:
                        logger.info(f"mod_type: {mod_type} identified as {matched_cv}")
                        formatted_mod_type_lst.append(matched_cv)
                        if matched_cv == "Delta":
                            mod_lv_dct[matched_cv] = self.raw_cv["Delta"].get("LEVEL", 0)
                            mass_shift_dct[matched_cv] = int(mod_type)
                            break
                        else:
                            mod_lv_dct[matched_cv] = max(
                                mod_lv_dct.get(matched_cv, 0), self.alias2cv[alia].get("LEVEL", 0)
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
        elif len(formatted_mod_type_lst) == len(mod_site_info_lst):
            mod_site_lst = [""] * len(formatted_mod_type_lst)
            formatted_mod_lst = zip(
                formatted_mod_type_lst, mod_site_lst, mod_site_info_lst
            )
        else:
            formatted_mod_lst = []

        mod_info_dct = {}
        if formatted_mod_lst:
            for mod_tp in formatted_mod_lst:
                mod_type = mod_tp[0]
                mod_order = self.raw_cv[mod_type].get("ORDER", 0)
                existed_mod_count = mod_info_dct.get(f'{mod_order}_{mod_type}', {}).get("MOD_COUNT", 0)
                existed_mod_site_lst = mod_info_dct.get(f'{mod_order}_{mod_type}', {}).get(
                    "MOD_SITE", []
                )
                existed_mod_site_info_lst = mod_info_dct.get(f'{mod_order}_{mod_type}', {}).get(
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
                elif mod_type == "Delta":
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
                else:
                    pass
                updated_mod_info = {
                    "MOD_CV": mod_type,
                    "MOD_LEVEL": mod_level,
                    "MOD_COUNT": mod_count,
                    "MOD_SITE": existed_mod_site_lst,
                    "MOD_SITE_INFO": existed_mod_site_info_lst,
                    "MOD_ORDER": mod_order
                }
                if mod_type in self.raw_cv:
                    updated_mod_info["MOD_ELEMENTS"] = self.raw_cv[mod_type].get("ELEMENTS", {})
                    if mod_type not in ["Delta", "DB"] and mod_type not in mass_shift_dct:
                        updated_mod_info["MOD_SHIFT"] = self.to_mass_shift(updated_mod_info["MOD_ELEMENTS"])
                    elif mod_type == "Delta" and mod_type in mass_shift_dct:
                        updated_mod_info["MOD_SHIFT"] = mass_shift_dct.get("Delta", 0)
                    else:
                        updated_mod_info["MOD_SHIFT"] = 0
                else:
                    raise ValueError(f"Unsupported modification type: {mod_type}")
                mod_info_dct[f'{mod_order}_{mod_type}'] = updated_mod_info
        else:
            pass
        mod_seg_levels_lst = []
        if mod_info_dct:
            for mod_seg in mod_info_dct:
                mod_seg_levels_lst.append(mod_info_dct[mod_seg].get("MOD_LEVEL", '0'))
        if mod_seg_levels_lst:
            max_mod_level = max(mod_seg_levels_lst)
            mod_seg_levels_str_lst = [str(i) for i in mod_seg_levels_lst]
            for mod_lv_str in mod_seg_levels_str_lst:
                if mod_lv_str.endswith(".1"):
                    max_mod_level += 0.1
                elif mod_lv_str.endswith(".2"):
                    max_mod_level += 0.2
                else:
                    pass
        else:
            max_mod_level = 0
        # mod_info_dct["MOD_LEVEL"] = max_mod_level

        return {"MOD_LEVEL": max_mod_level, "MOD_INFO": mod_info_dct}

    def format_residue(self, info: dict) -> dict:
        residue_info_dct = {}
        link_lst = info.get("LINK", [""])
        if link_lst:
            link = link_lst[0]
        else:
            link = ""
        num_o_lst = info.get("NUM_O", ["0"])
        num_o = 0
        if num_o_lst:
            if num_o_lst[0] != "":
                num_o_str = str(num_o_lst[0]).strip("O")
                num_o = int(num_o_str)
            else:
                pass
        else:
            pass

        residue_info_dct["LINK"] = link
        residue_info_dct["MOD"] = self.format_mods(info)
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
