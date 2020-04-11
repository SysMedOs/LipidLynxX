# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
from typing import Dict, List, Union

from jsonschema import Draft7Validator, RefResolver
from natsort import natsorted
import regex as re

from lynx.utils.params_loader import load_output_rule
from lynx.models.defaults import (
    api_version,
    lynx_schema_cfg,
    core_schema,
    core_schema_path,
    mod_db_level_lst,
    default_output_rules,
)
from lynx.utils.file_readers import get_abs_path
from lynx.utils.log import logger
from lynx.utils.toolbox import check_json


class Modifications(object):
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
        self.schema = schema
        self.type = "Modification"
        self.mod_level = str(mod_info.get("MOD_LEVEL", 0))
        with open(get_abs_path(lynx_schema_cfg[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(
                json.load(s_obj),
                resolver=RefResolver(
                    f"file://{core_schema_path}", referrer=core_schema
                ),
            )

        self.db_count = db
        self.sum_mod_info = self.to_sum_info()

        self.mod_id = self.sum_mod_info.get("id", "")
        self.mod_linked_ids = self.sum_mod_info.get("linked_ids", {})
        self.mod_list = self.sum_mod_info.get("info", {})

    def __str__(self):
        return self.to_json()

    def __repr__(self):
        return self.to_json()

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
            if self.mod_info[mod].get("MOD_CV", "") not in ["", "DB"]:
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
                if o in mod_seg_lst or o in self.mod_separators:
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
                    elif o.upper().endswith("_SEPARATOR"):
                        mod_seg_str += self.mod_separators.get(o, "")
                    elif re.search("BRACKET", o.upper()):
                        mod_seg_str += self.mod_separators.get(o, "")
                    else:
                        mod_seg_str += str(mod_dct.get(o, ""))

            mod_seg_str = re.sub(r"\{\}", "", mod_seg_str)
            mod_seg_str = re.sub(r",,", ",", mod_seg_str)
            mod_str_dct[mod] = mod_seg_str.strip(",")
        if get_db_only and db_idx and mod_str_dct:
            mod_str_dct = {"0": mod_str_dct[db_idx]}
        else:
            if db_idx in mod_str_dct:
                del mod_str_dct[db_idx]
            else:
                pass
        mod_str_lst = []
        if mod_str_dct:
            mod_str_dct_idx = natsorted(list(mod_str_dct.keys()))
            for k in mod_str_dct_idx:
                mod_str_lst.append(mod_str_dct[k])
        return ",".join(mod_str_lst).strip(",")

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
        if float(level) > float(self.mod_level):
            raise ValueError(
                f'Cannot convert to higher level than the mod_level "{self.mod_level}". Input:{level}'
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

        # add DB part
        if level.endswith(".1"):
            db_str = self.to_mod_site(get_db_only=True)
            if db_str:
                mod_str = ",".join([db_str, mod_str]).strip(",")
            else:
                pass
        elif level.endswith(".2"):
            db_str = self.to_mod_site_info(get_db_only=True)
            if db_str:
                mod_str = ",".join([db_str, mod_str]).strip(",")
            else:
                pass
        else:
            pass
        mod_str = mod_str.strip(",")

        return mod_str

    def to_all_levels(self, as_list: bool = False) -> Union[Dict[str, str], List[str]]:
        all_levels_dct = {}
        if self.mod_level in mod_db_level_lst:
            mod_idx = mod_db_level_lst.index(self.mod_level)
            output_levels_lst = mod_db_level_lst[: mod_idx + 1]
        else:
            raise ValueError(f"Modification level not supported: {self.mod_level}")
        if self.mod_level.endswith(".2"):
            pass
        elif self.mod_level.endswith(".1"):
            output_levels_lst = [
                o_lv for o_lv in output_levels_lst if not o_lv.endswith(".2")
            ]
        else:
            output_levels_lst = [
                o_lv
                for o_lv in output_levels_lst
                if not o_lv.endswith(".1") and not o_lv.endswith(".2")
            ]
        for level in output_levels_lst:
            all_levels_dct[level] = self.to_mod_level(level)
        all_levels_info = all_levels_dct

        return all_levels_info

    def get_mod_info(self) -> list:
        mod_js_lst = []
        for mod_idx in self.mod_info:
            mod_seg_dct = self.mod_info[mod_idx]
            mod_sites_lst = mod_seg_dct.get("MOD_SITE", [])
            if mod_sites_lst:
                try:
                    sites = [int(i) for i in mod_seg_dct.get("MOD_SITE", [])]
                except ValueError:
                    sites = []
            else:
                sites = []
            mod_js_dct = {
                "cv": mod_seg_dct.get("MOD_CV", ""),
                "count": mod_seg_dct.get("MOD_COUNT", 0),
                "site": sites,
                "site_info": mod_seg_dct.get("MOD_SITE_INFO", []),
                "order": mod_seg_dct.get("MOD_ORDER", 0),
                "elements": mod_seg_dct.get("MOD_ELEMENTS", {}),
                "mass_shift": mod_seg_dct.get("MOD_MASS_SHIFT", 0),
            }
            mod_js_lst.append(mod_js_dct)
        return mod_js_lst

    def to_sum_info(self):
        linked_ids = self.to_all_levels()
        mod_id = linked_ids.get(self.mod_level, "")
        if float(self.mod_level) > 0 and mod_id:
            sum_mod_info_dct = {
                "api_version": api_version,
                "type": self.type,
                "id": mod_id,
                "level": self.mod_level,
                "linked_ids": linked_ids,
                "linked_levels": natsorted(list(linked_ids.keys())),
                "info": self.get_mod_info(),
            }
        elif float(self.mod_level) == 0:
            sum_mod_info_dct = {}
        else:
            raise ValueError(
                f"Cannot format_mods modification code to level {self.mod_level} "
                f"from input: {self.mod_info}"
            )

        return sum_mod_info_dct

    def to_json(self):
        mod_json_str = json.dumps(self.sum_mod_info)

        if check_json(validator=self.validator, json_obj=json.loads(mod_json_str)):
            return mod_json_str
        else:
            raise Exception(f"Schema test FAILED. Schema {self.schema}")


def merge_mods(
    mods_collection: Union[dict, list],
    db: int = 0,
    schema: str = "lynx_mod",
    output_rules: dict = default_output_rules,
    nomenclature: str = "LipidLynxX",
) -> Modifications:
    sum_mods_dct = {}
    if isinstance(mods_collection, list):
        pass
    elif isinstance(mods_collection, dict):
        mods_collection = [
            mods_collection[mod] for mod in mods_collection
        ]  # type: list
    else:
        raise TypeError(
            f"Requires multiple Mods information in list or dict, "
            f"got type: {type(mods_collection)} for {mods_collection}"
        )

    for mod_seg in mods_collection:
        if isinstance(mod_seg, Modifications):
            mod_info = mod_seg.mod_info
        elif isinstance(mod_seg, dict):
            mod_info = mod_seg.get("MOD_INFO", {})
        else:
            raise TypeError(
                f"Requires multiple Mods information in dict or as Mods object, "
                f"got type: {type(mod_seg)} for {mod_seg}"
            )

        for mod_idx in mod_info:
            if re.search(r"DB", mod_idx):
                pass
            else:
                if mod_idx not in sum_mods_dct:
                    sum_mods_dct[mod_idx] = mod_info[mod_idx]
                else:
                    existed_count = sum_mods_dct[mod_idx].get("MOD_COUNT", 0)
                    sum_mods_dct[mod_idx]["MOD_COUNT"] = (
                        mod_info[mod_idx].get("MOD_COUNT", 0) + existed_count
                    )
    max_level = 0
    for sum_mod_idx in sum_mods_dct:
        sum_mod_seg_info = sum_mods_dct[sum_mod_idx]
        mod_seg_count = sum_mod_seg_info["MOD_COUNT"]
        if sum_mod_seg_info["MOD_LEVEL"] > 3:
            sum_mod_seg_info["MOD_LEVEL"] = 3
            max_level = max(max_level, sum_mod_seg_info["MOD_LEVEL"])
        if sum_mod_seg_info["MOD_SITE"]:
            sum_mod_seg_info["MOD_SITE"] = [""] * mod_seg_count
        if sum_mod_seg_info["MOD_SITE_INFO"]:
            sum_mod_seg_info["MOD_SITE_INFO"] = [""] * mod_seg_count
        sum_mods_dct[sum_mod_idx] = sum_mod_seg_info

    sum_mod_obj = Modifications({"MOD_LEVEL": max_level, "MOD_INFO": sum_mods_dct})

    return sum_mod_obj


if __name__ == "__main__":

    usr_mod_info = {
        "MOD_LEVEL": 5.2,
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
                "MOD_COUNT": 1,
                "MOD_SITE": ["12"],
                "MOD_SITE_INFO": ["R"],
                "MOD_ORDER": 5.01,
                "MOD_ELEMENTS": {"O": 1},
                "MOD_MASS_SHIFT": 16,
            },
        },
    }

    usr_mod_obj = Modifications(usr_mod_info)
    logger.debug(usr_mod_obj)
    mod_json = usr_mod_obj.to_json()

    usr_sum_mods_obj = merge_mods([usr_mod_obj, usr_mod_obj])

    logger.debug(usr_sum_mods_obj)
    sum_mod_json = usr_sum_mods_obj.to_json()

    logger.info("FINISHED")
