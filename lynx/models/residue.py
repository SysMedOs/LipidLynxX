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

import json
import os

from jsonschema import Draft7Validator, RefResolver
import regex as re

from lynx.utils.params_loader import load_output_rule
from lynx.models.defaults import res_schema, res_schema_path, default_output_rules
from lynx.models.modification import Modifications, merge_mods
from lynx.utils.log import app_logger
from lynx.utils.toolbox import check_json


class Residue(object):
    def __init__(
        self,
        residue_info: dict,
        schema: str = "lynx_residues",
        output_rules: dict = default_output_rules,
        nomenclature: str = "LipidLynxX",
        logger=app_logger,
    ):
        self.logger = logger
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

        self.mod_obj = Modifications(mod_info, nomenclature=nomenclature)
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
                            o_seg_rgx = self.res_rule.get("RESIDUE", {}).get("NUM_O")
                            if o_seg_rgx:
                                if num_o == 1:
                                    if re.match(o_seg_rgx, str(num_o)):
                                        res_str += str(num_o)
                                    elif re.match(o_seg_rgx, "1"):
                                        res_str += str("1")
                                    elif re.match(o_seg_rgx, "O"):
                                        res_str += str("O")
                                    else:
                                        res_str += str(num_o)
                                else:
                                    if re.match(o_seg_rgx, str(num_o)):
                                        res_str += str(num_o)
                                    elif re.match(o_seg_rgx, f"{num_o}O"):
                                        res_str += f"{num_o}O"
                                    elif re.match(o_seg_rgx, f"O{num_o}"):
                                        res_str += f"O{num_o}"
                                    else:
                                        res_str += str(num_o)
                            else:
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
        if check_json(self.validator, json.loads(fa_json_str, logger=self.logger)):
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
        all_residues[rm].get("MOD", {"MOD_LEVEL": 0, "MOD_INFO": {}})
        for rm in all_residues
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

    sum_res_dct["MOD"] = {
        "MOD_LEVEL": sum_mods_obj.mod_level,
        "MOD_INFO": sum_mods_obj.mod_info,
    }

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
    app_logger.debug(usr_res_obj.linked_ids)
    # usr_res_json = res_obj.to_json()
    app_logger.info("FINISHED")
