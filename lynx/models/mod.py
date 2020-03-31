# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os
from typing import Dict, List, Union

from jsonschema import Draft7Validator, RefResolver
from natsort import natsorted

from lynx.controllers.params_loader import load_output_rule
from lynx.models.cv import CV
from lynx.models.defaults import (
    api_version,
    cv_elements_info,
    elem_nominal_info,
    lynx_schema_cfg,
    core_schema,
    core_schema_path,
    mod_db_level_lst,
    default_cv_file,
    default_output_rules,
)
from lynx.utils.file_readers import get_abs_path
from lynx.utils.log import logger
from lynx.utils.toolbox import check_json


class Mods(object):
    def __init__(
        self,
        mod_info: dict,
        db: int = 0,
        schema: str = "lynx_mod",
        output_rules: dict = default_output_rules,
        nomenclature: str = "LipidLynxX",
    ):
        self.export_rule = load_output_rule(output_rules, nomenclature)
        self.mod_rule = self.export_rule.get("MODS", None)
        self.mod_rule_orders = self.mod_rule.get("MOD", {}).get("ORDER", [])
        self.mod_separators = self.export_rule.get("SEPARATORS", [])
        if not self.mod_rule:
            raise ValueError(
                f"Cannot find output rule for 'MODS' from nomenclature: {nomenclature}."
            )
        self.mod_info = mod_info.get("MOD_INFO", {})
        self.schema = "lynx_mod"
        self.type = "Modification"
        self.max_mod_level = str(mod_info.get("MOD_LEVEL", 0))
        with open(get_abs_path(lynx_schema_cfg[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(
                json.load(s_obj),
                resolver=RefResolver(
                    f"file://{core_schema_path}", referrer=core_schema
                ),
            )

        self.db_count = db
        self.sum_mod_info = self.get_sum_info()

        self.mod_id = self.sum_mod_info.get("id", "")
        self.mod_linked_ids = self.sum_mod_info.get("linked_ids", {})
        self.mod_list = self.sum_mod_info.get("info", {})

        # logger.info(
        #     f"Level {self.mod_level:3s} modification created from: {self.mod_code}"
        # )

    @staticmethod
    def __add_angle_brackets__(mod_str: str) -> str:
        if mod_str and mod_str not in ["+0", "-0", "0"]:
            return f"<{mod_str}>"
        else:
            return ""

    def get_sum_elements_shift(self):
        sum_mod_elem_dct = {}
        for mod in self.mod_info:
            mod_cv = mod["MOD_CV"]
            mod_count = mod["MOD_COUNT"]
            mod_elem_dct = mod["MOD_ELEMENTS"]
            if mod_cv != "DB":
                for elem in mod_elem_dct:
                    if elem in sum_mod_elem_dct:
                        sum_mod_elem_dct[elem] += mod_count * mod_elem_dct[elem]
                    else:
                        sum_mod_elem_dct[elem] = mod_count * mod_elem_dct[elem]

        return sum_mod_elem_dct

    def get_sum_info(self):
        linked_ids = self.to_all_levels()
        mod_id = linked_ids.get(self.max_mod_level, "")
        if mod_id:
            if mod_id == self.mod_code:
                sum_mod_info_dct = {
                    "api_version": api_version,
                    "type": self.type,
                    "id": self.mod_code,
                    "level": self.max_mod_level,
                    "all_levels": self.max_mod_level,
                    "linked_ids": self.to_all_levels(),
                    "info": self.mod_info,
                }
            else:
                sum_mod_info_dct = {
                    "api_version": api_version,
                    "type": self.type,
                    "id": mod_id,
                    # "input_name": self.mod_code,
                    "level": self.max_mod_level,
                    "linked_ids": self.to_all_levels(),
                    "info": self.mod_info,
                }
        else:
            raise ValueError(
                f"Cannot format_mods modification code to level {self.max_mod_level} "
                f"from input: {self.mod_code}"
            )

        return sum_mod_info_dct

    def to_json(self):
        mod_json_str = json.dumps(self.sum_mod_info)

        if check_json(validator=self.validator, json_obj=json.loads(mod_json_str)):
            return mod_json_str
        else:
            raise Exception(f"Schema test FAILED. Schema {self.schema}")

    def to_mass_shift(self) -> str:
        mass_shift = 0
        for mod in self.mod_info:
            mass_shift += self.mod_info[mod].get("MOD_MASS_SHIFT", 0)

        return f"{mass_shift:+}"

    def to_elements(self):
        sum_elements = {}
        mod_elem_lst = ["C", "O", "N", "S", "H", "Na"]
        mod_str_lst = []
        for mod in self.mod_info:
            mod_elements = self.mod_info[mod].get("MOD_ELEMENTS", 0)
            for elem in mod_elements:
                sum_elements[elem] = sum_elements.get(elem, 0) + mod_elements.get(
                    elem, 0
                )

        for mod_elem in mod_elem_lst:
            if mod_elem in sum_elements:
                mod_elem_count = f"{sum_elements[mod_elem]:+}"
                if mod_elem_count == "+1":
                    mod_elem_count = "+"
                elif mod_elem_count == "-1":
                    mod_elem_count = "-"
                mod_str_lst.append(f"{mod_elem_count}{mod_elem}")

        return f'{",".join(mod_str_lst)}'

    def to_mod_base(
        self,
        mod_seg_lst: Union[list, tuple] = ("MOD_COUNT", "MOD_CV"),
        get_db_only: bool = False,
    ) -> str:
        mod_str_dct = {}
        db_idx = None
        for mod in self.mod_info:
            mod_seg_str = ""
            is_mod_sites_with_info = False
            mod_sites_lst = []
            mod_dct = self.mod_info[mod]
            cv = mod_dct.get("MOD_CV", None)
            if cv in ["", "DB"]:
                db_idx = mod
            for o in self.mod_rule_orders:
                if o in mod_seg_lst:
                    if o == "MOD_COUNT":
                        mod_count = mod_dct.get(o, 1)
                        if mod_count > 1 and cv not in ["", "DB"]:
                            mod_seg_str += str(mod_count)
                        elif mod_count == 1:
                            pass
                        elif cv in ["", "DB"]:
                            pass
                        else:
                            raise ValueError(
                                f"Modification count must >= 1, got value: {mod_count}"
                            )
                    elif o == "MOD_CV":
                        if cv in ["", "DB"]:
                            pass
                        else:
                            mod_seg_str += cv
                    elif o == "MOD_SITE" and "MOD_SITE_INFO" not in mod_seg_lst:
                        mod_sites_lst = mod_dct.get("MOD_SITE", [])
                        if mod_sites_lst:
                            mod_seg_str += ",".join(mod_sites_lst)
                        else:
                            pass
                    elif o == "MOD_SITE" and "MOD_SITE_INFO" in mod_seg_lst:
                        mod_sites_lst = mod_dct.get("MOD_SITE", [])
                        is_mod_sites_with_info = True
                    elif o == "MOD_SITE_INFO" and is_mod_sites_with_info:
                        mod_sites_info_lst = mod_dct.get("MOD_SITE_INFO", [])
                        mod_sites_cmb_lst = zip(mod_sites_lst, mod_sites_info_lst)
                        mod_sites_str_lst = [f"{t[0]}{t[1]}" for t in mod_sites_cmb_lst]
                        if mod_sites_str_lst:
                            mod_seg_str += ",".join(mod_sites_str_lst)
                        else:
                            pass
                    elif o == "SITE_BRACKET_LEFT":
                        mod_seg_str += self.mod_separators.get("SITE_BRACKET_LEFT", "{")
                    elif o == "SITE_BRACKET_RIGHT":
                        mod_seg_str += self.mod_separators.get("SITE_BRACKET_RIGHT", "}")
                    else:
                        mod_seg_str += str(mod_dct.get(o, ""))
            mod_str_dct[mod] = mod_seg_str
        if get_db_only and db_idx and mod_str_dct:
            mod_str_dct = {"0": mod_str_dct[db_idx]}
        mod_str_lst = []
        if mod_str_dct:
            mod_str_dct_idx = natsorted(list(mod_str_dct.keys()))
            for k in mod_str_dct_idx:
                mod_str_lst.append(mod_str_dct[k])
        return ",".join(mod_str_lst).strip(',')

    def to_mod_count(self, get_db_only: bool = False):
        return self.to_mod_base(
            mod_seg_lst=["MOD_COUNT", "MOD_CV"], get_db_only=get_db_only
        )

    def to_mod_site(self, get_db_only: bool = False):
        return self.to_mod_base(
            mod_seg_lst=[
                "MOD_COUNT",
                "MOD_CV",
                "SITE_BRACKET_LEFT",
                "MOD_SITE",
                "SITE_BRACKET_RIGHT",
            ],
            get_db_only=get_db_only,
        )

    def to_mod_site_info(self, get_db_only: bool = False):
        return self.to_mod_base(
            mod_seg_lst=[
                "MOD_COUNT",
                "MOD_CV",
                "SITE_BRACKET_LEFT",
                "MOD_SITE",
                "MOD_SITE_INFO",
                "SITE_BRACKET_RIGHT",
            ],
            get_db_only=get_db_only,
        )

    def to_mod_level(self, level: Union[int, float, str] = 0) -> str:
        mod_str = ""
        if not isinstance(level, str):
            level = str(level)
        if float(level) > float(self.max_mod_level):
            raise ValueError(
                f'Cannot convert to higher level than the mod_level "{self.max_mod_level}". Input:{level}'
            )

        if level.startswith("0"):
            mod_str = ""
        elif level.startswith("1"):
            mod_str = self.to_mass_shift()
        elif level.startswith("2"):
            mod_str = self.to_elements()
        elif level.startswith("3"):
            mod_str = self.to_mod_count()
        elif level.startswith("4"):
            mod_str = self.to_mod_site()
        elif level.startswith("5"):
            mod_str = self.to_mod_site_info()
        else:
            raise ValueError(f"Currently not supported modification level: {level}")
        if float(level) < 4:
            if level.endswith(".1"):
                db_str = self.to_mod_site(get_db_only=True)
                mod_str = ','.join([db_str, mod_str])
            elif level.endswith(".2"):
                db_str = self.to_mod_site_info(get_db_only=True)
                mod_str = ','.join([db_str, mod_str])
            else:
                pass
        else:
            pass
        mod_str = mod_str.strip(",")

        return mod_str

    def to_all_levels(self, as_list: bool = False) -> Union[Dict[str, str], List[str]]:
        all_levels_dct = {}
        if self.max_mod_level in mod_db_level_lst:
            mod_idx = mod_db_level_lst.index(self.max_mod_level)
            output_levels_lst = mod_db_level_lst[: mod_idx + 1]
        else:
            raise ValueError(f"Modification level not supported: {self.max_mod_level}")

        for level in output_levels_lst:
            all_levels_dct[level] = self.to_mod_level(level)
        all_levels_info = all_levels_dct

        return all_levels_info


if __name__ == "__main__":
    usr_mod_info = {
        "MOD_LEVEL": 3.2,
        "MOD_INFO": {
            "0.01_DB": {
                "MOD_CV": "DB",
                "MOD_LEVEL": 0.2,
                "MOD_COUNT": 2,
                "MOD_SITE": ["9", "11"],
                "MOD_SITE_INFO": ["Z", "Z"],
                "MOD_ORDER": 0.01,
                "MOD_ELEMENTS": {"H": -2},
                "MOD_MASS_SHIFT": 0,
            },
            "5.01_OH": {
                "MOD_CV": "OH",
                "MOD_LEVEL": 5,
                "MOD_COUNT": 2,
                "MOD_SITE": [],
                "MOD_SITE_INFO": [],
                "MOD_ORDER": 5.01,
                "MOD_ELEMENTS": {"O": 1},
                "MOD_MASS_SHIFT": 16,
            },
        },
    }
    # usr_mod_info = {
    #     "MOD_LEVEL": 5.2,
    #     "MOD_INFO": {
    #         "0.01_DB": {
    #             "MOD_CV": "DB",
    #             "MOD_LEVEL": 0.2,
    #             "MOD_COUNT": 2,
    #             "MOD_SITE": ["9", "11"],
    #             "MOD_SITE_INFO": ["Z", "Z"],
    #             "MOD_ORDER": 0.01,
    #             "MOD_ELEMENTS": {"H": -2},
    #             "MOD_MASS_SHIFT": 0,
    #         },
    #         "5.01_OH": {
    #             "MOD_CV": "OH",
    #             "MOD_LEVEL": 5,
    #             "MOD_COUNT": 1,
    #             "MOD_SITE": ["12"],
    #             "MOD_SITE_INFO": ["R"],
    #             "MOD_ORDER": 5.01,
    #             "MOD_ELEMENTS": {"O": 1},
    #             "MOD_MASS_SHIFT": 16,
    #         },
    #     },
    # }

    mod_obj = Mods(usr_mod_info)
    mod_obj

    logger.info("FINISHED")
