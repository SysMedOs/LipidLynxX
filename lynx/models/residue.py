# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os

from jsonschema import Draft7Validator, RefResolver
import regex as re

from lynx.utils.params_loader import load_output_rule
from lynx.models.defaults import (
    api_version,
    mod_level_lst,
    hg_schema,
    hg_schema_path,
    res_schema,
    res_schema_path,
    default_output_rules,
)
from lynx.models.modifications import Modifications
from lynx.models.mod import Mods, merge_mods
from lynx.models.patterns import fa_rgx
from lynx.utils.log import logger
from lynx.utils.toolbox import check_json


class LipidClass(object):
    def __init__(self, hg_code: str):
        self.hg = hg_code
        self.schema = "lynx_hg"
        self.lynx_type = "HeadGroup"
        resolver = RefResolver(
            referrer=hg_schema, base_uri=f"file://{os.path.dirname(hg_schema_path)}/"
        )
        self.validator = Draft7Validator(hg_schema, resolver=resolver)

        self.sum_info = {
            "api_version": api_version,
            "id": hg_code,
            "type": self.lynx_type,
            "level": "S",
        }
        self.id = hg_code

    def to_json(self):
        mod_json_str = json.dumps(self.sum_info)
        if check_json(validator=self.validator, json_obj=json.loads(mod_json_str)):
            return mod_json_str
        else:
            raise Exception(f"Schema test FAILED. Schema {self.schema}")


class FattyAcid(object):
    def __init__(self, lipid_code: str, db: int = 0):

        self.lipid_code = lipid_code.strip("FA")
        self.schema = "lynx_fa"
        self.type = "FattyAcid"
        resolver = RefResolver(
            referrer=res_schema, base_uri=f"file://{os.path.dirname(res_schema_path)}/"
        )
        self.validator = Draft7Validator(res_schema, resolver=resolver)

        self.fa_info_dct = self.__post_init__()
        self.fa_info_dct["id"] = self.lipid_code
        self.info = self.fa_info_dct["info"]

        self.mod_info = self.fa_info_dct.get("mod_obj", None)
        if db and db < self.fa_info_dct["info"]["db"]:
            self.db_count = db
        else:
            self.db_count = self.fa_info_dct["info"]["db"]
        self.id = self.fa_info_dct.get("id", "")
        self.fa_level = self.fa_info_dct.get("level", "")
        self.is_modified = self.fa_info_dct["info"].get("is_modified", False)
        self.fa_linked_ids = self.fa_info_dct.get("linked_ids", {})
        logger.debug(
            f"Level {self.fa_level:4s} FattyAcid created from: {self.lipid_code}"
        )

    def __post_init__(self):

        fa_info_dct = {"api_version": api_version, "type": self.type, "info": {}}
        fa_match = fa_rgx.match(self.lipid_code)

        is_modified = False

        if fa_match:

            fa_matched_dct = fa_match.groupdict()
            mod_code = fa_matched_dct.get("mod", None)
            fa_link = fa_matched_dct.get("link", "")
            if not fa_link:
                fa_link = "FA"
            fa_info_dct["info"] = {
                "c": int(fa_matched_dct.get("c", 0)),
                "db": int(fa_matched_dct.get("db", 0)),
                "link": fa_link,
            }
            fa_seg_str = (
                f'{fa_link}{fa_matched_dct.get("c", 0)}:{fa_matched_dct.get("db", 0)}'
            )
            if fa_seg_str.lower().startswith("none"):
                fa_seg_str = fa_seg_str[4:]
            if mod_code and mod_code.strip("<>"):  # mod_code can be None or '<>'
                mod_obj = Modifications(mod_code)
                fa_info_dct["mod_obj"] = mod_obj
                fa_info_dct["level"] = f"{mod_obj.mod_level}"
                fa_linked_ids_dct = {}
                mod_linked_ids_dct = mod_obj.mod_linked_ids  # type: dict
                for lv in mod_linked_ids_dct:
                    fa_linked_ids_dct[f"{lv}"] = f"{fa_seg_str}{mod_linked_ids_dct[lv]}"
                fa_info_dct["linked_ids"] = fa_linked_ids_dct

                if mod_obj.mod_list:
                    for mod_seg in mod_obj.mod_list:
                        if (
                            mod_seg.get("cv", "DB") != "DB"
                            and mod_seg.get("count", 0) != 0
                        ):
                            is_modified = True
                            fa_info_dct["info"]["is_modified"] = is_modified
                            fa_info_dct["info"]["modifications"] = mod_obj.sum_mod_info
                        else:
                            fa_info_dct["info"]["is_modified"] = is_modified
            else:
                if fa_matched_dct.get("db", 0) == 0:
                    fa_info_dct["level"] = "5.2"
                    for mod_lv in mod_level_lst:
                        fa_info_dct["linked_ids"][mod_lv] = fa_seg_str
                else:
                    fa_info_dct["level"] = "5"
                    fa_info_dct["linked_ids"] = {
                        "0": fa_seg_str,
                        "1": fa_seg_str,
                        "2": fa_seg_str,
                        "3": fa_seg_str,
                        "4": fa_seg_str,
                        "5": fa_seg_str,
                    }
                if self.lipid_code != fa_seg_str:
                    self.lipid_code = fa_seg_str
                fa_info_dct["info"]["is_modified"] = is_modified
            return fa_info_dct
        else:
            raise ValueError(f"Cannot parse FA sting: {self.lipid_code}")

    def to_json(self):
        fa_lite_info_dct = self.fa_info_dct
        fa_lite_info_dct.pop("mod_obj", None)
        fa_json_str = json.dumps(fa_lite_info_dct)
        if check_json(self.validator, json.loads(fa_json_str)):
            return fa_json_str
        else:
            raise Exception(f"JSON Schema check FAILED. Schema {self.schema}")

    def to_segments(self, mod_level: str):

        out_fa_info_dct = self.info
        if (
            self.mod_info
            and self.mod_info.mod_linked_ids
            and mod_level in mod_level_lst
            and float(mod_level) <= float(mod_level)
        ):
            out_fa_info_dct["mod_id"] = self.mod_info.mod_linked_ids.get(mod_level, "")
        else:
            out_fa_info_dct["mod_id"] = ""
        return out_fa_info_dct


class Residue(object):
    def __init__(
        self,
        residue_info: dict,
        schema: str = "lynx_residues",
        output_rules: dict = default_output_rules,
        nomenclature: str = "LipidLynxX",
    ):
        self.export_rule = load_output_rule(output_rules, nomenclature)
        self.res_rule = self.export_rule.get("RESIDUES", None)
        self.res_rule_orders = self.res_rule.get("RESIDUE", {}).get("ORDER", [])
        self.res_separators = self.export_rule.get("SEPARATORS", [])
        self.res_info = residue_info
        self.__replace_mdt__()

        self.schema = schema
        self.type = "FattyAcid"
        resolver = RefResolver(
            referrer=res_schema, base_uri=f"file://{os.path.dirname(res_schema_path)}/"
        )
        self.validator = Draft7Validator(res_schema, resolver=resolver)

        mod_info = residue_info.get("MOD", {})

        self.mod_obj = Mods(mod_info)
        self.sum_mod_info = self.mod_obj.sum_mod_info
        self.mod_level = self.mod_obj.mod_level
        self.res_level = self.mod_level
        if float(self.mod_level) > 0:
            self.is_modified = True
        else:
            self.is_modified = False

        if self.is_modified and self.sum_mod_info:
            self.linked_levels = self.sum_mod_info.get("linked_levels", ["0"])
        else:
            self.linked_levels = ["0"]

        self.linked_ids = self.__post_init__()

    def __replace_mdt__(self):
        link = self.res_info.get("LINK", "")
        if link.lower() == "m":
            self.res_info["LINK"] = ""
            self.res_info["NUM_O"] = 1 + self.res_info.get("NUM_O", 0)
        elif link.lower() == "d":
            self.res_info["LINK"] = ""
            self.res_info["NUM_O"] = 2 + self.res_info.get("NUM_O", 0)
        elif link.lower() == "t":
            self.res_info["LINK"] = ""
            self.res_info["NUM_O"] = 3 + self.res_info.get("NUM_O", 0)
        else:
            pass

    def __post_init__(self):
        res_str_dct = {}
        num_o = self.res_info.get("NUM_O", 0)
        for lv in self.linked_levels:
            res_str = ""
            for o in self.res_rule_orders:
                if o in self.res_info or o in self.res_separators or o in ["SUM_MODS"]:
                    if o in ["SUM_MODS", "MODS", "MOD"]:
                        res_str += self.sum_mod_info.get("linked_ids", {}).get(lv, "")
                    elif o == "NUM_O":
                        if num_o > 0:
                            res_str += str(num_o)
                        else:
                            pass
                    elif o.upper().endswith("_SEPARATOR"):
                        res_str += self.res_separators.get(o, "")
                        if num_o == 0:
                            res_str = res_str.strip(
                                self.res_separators.get("O_SEPARATOR", "")
                            )
                    elif re.search("BRACKET", o.upper()):
                        res_str += self.res_separators.get(o, "")
                    else:
                        res_str += str(self.res_info.get(o, ""))
            na_brackets_lst = [r"\<\>", r"\{\}", r"\[\]", r"\(\)"]
            for b in na_brackets_lst:
                res_str = re.sub(b, "", res_str)
            res_str_dct[lv] = res_str.strip(";")

        return res_str_dct

    def to_json(self):
        fa_lite_info_dct = self.fa_info_dct
        fa_lite_info_dct.pop("mod_obj", None)
        fa_json_str = json.dumps(fa_lite_info_dct)
        if check_json(self.validator, json.loads(fa_json_str)):
            return fa_json_str
        else:
            raise Exception(f"JSON Schema check FAILED. Schema {self.schema}")


def merge_residues(
    all_residues: dict,
    schema: str = "lynx_residues",
    output_rules: dict = default_output_rules,
    nomenclature: str = "LipidLynxX",
) -> Residue:

    sum_res_dct = {}
    if isinstance(all_residues, dict):
        pass
    else:
        raise TypeError(
            f"Requires multiple Residues in dict, "
            f"got type: {type(all_residues)} for {all_residues}"
        )

    all_mod_lst = [
        all_residues[rm].get("MOD", {"MOD_LEVEL": 0, "MOD_INFO": {}}) for rm in all_residues
    ]
    sum_mods_obj = merge_mods(all_mod_lst)

    for res in all_residues:
        res_info = all_residues[res]
        for res_seg in res_info:
            if re.search(r"MOD", res_seg):
                pass
            else:
                if res_seg not in sum_res_dct:
                    sum_res_dct[res_seg] = res_info[res_seg]
                else:
                    existed_count = sum_res_dct.get(res_seg, None)
                    res_seg_count = res_info.get(res_seg, None)
                    if res_seg_count:
                        if isinstance(existed_count, int) and isinstance(
                            res_seg_count, int
                        ):
                            sum_res_dct[res_seg] = res_seg_count + existed_count
                        elif isinstance(existed_count, str) and isinstance(
                            res_seg_count, str
                        ):
                            sum_res_dct[res_seg] = res_seg_count + existed_count
                        else:
                            raise TypeError
                    else:
                        pass

    sum_res_dct["MOD"] = {"MOD_LEVEL": sum_mods_obj.mod_level, "MOD_INFO": sum_mods_obj.mod_info}

    sum_res_obj = Residue(sum_res_dct, schema, output_rules, nomenclature)

    return sum_res_obj


if __name__ == "__main__":

    usr_res_info = {
        "d18:1": {
            "LINK": "d",
            "MOD": {"MOD_LEVEL": 0, "MOD_INFO": {}},
            "NUM_C": 18,
            "NUM_DB": 1,
            "NUM_O": 0,
        },
        "18:2(9Z,11Z)(12OH)": {
            "LINK": "",
            "MOD": {
                "MOD_LEVEL": 4.2,
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
                        "MOD_LEVEL": 4,
                        "MOD_COUNT": 1,
                        "MOD_SITE": ["12"],
                        "MOD_SITE_INFO": [""],
                        "MOD_ORDER": 5.01,
                        "MOD_ELEMENTS": {"O": 1},
                        "MOD_MASS_SHIFT": 16,
                    },
                },
            },
            "NUM_C": 18,
            "NUM_DB": 2,
            "NUM_O": 0,
        },
    }
    # for r in usr_res_info:
    #     usr_res_obj = Residue(usr_res_info[r])
    #     logger.debug(usr_res_obj.linked_ids)
    #     # usr_res_json = res_obj.to_json()

    usr_res_obj = merge_residues(usr_res_info)
    logger.debug(usr_res_obj.linked_ids)
    # usr_res_json = res_obj.to_json()
    logger.info("FINISHED")
