# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json

from jsonschema import Draft7Validator

from lipidlynx.controllers.general_functions import check_json, get_abs_path

from lipidlynx.models.defaults import lynx_schema
from lipidlynx.models.log import logger
from lipidlynx.models.modification import Modifications
from lipidlynx.models.patterns import fa_rgx


class HeadGroup(object):
    def __init__(self, hg_code: str):
        self.hg = hg_code
        self.schema = "lynx_hg"
        with open(get_abs_path(lynx_schema[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(json.load(s_obj))

        self.sum_info = {"hg_id": hg_code}
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
        with open(get_abs_path(lynx_schema[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(json.load(s_obj))

        self.fa_info_dct = self.__post_init__()
        self.fa_info_dct["fa_id"] = self.lipid_code

        self.mod_info = self.fa_info_dct.get("mod_obj", None)
        if db and db < self.fa_info_dct["fa_info"]["db"]:
            self.db_count = db
        else:
            self.db_count = self.fa_info_dct["fa_info"]["db"]
        self.id = self.fa_info_dct.get("fa_id", "")
        self.fa_level = self.fa_info_dct.get("fa_level", "")
        self.is_modified = self.fa_info_dct["fa_info"].get("is_modified", False)
        self.fa_linked_ids = self.fa_info_dct.get("fa_linked_ids", {})
        logger.info(
            f"Level {self.fa_level:4s} FattyAcid created from: {self.lipid_code}"
        )

    def __post_init__(self):

        fa_info_dct = {}
        fa_match = fa_rgx.match(self.lipid_code)

        is_modified = False

        if fa_match:
            fa_matched_dct = fa_match.groupdict()
            mod_code = fa_matched_dct.get("mod", None)
            fa_link = fa_matched_dct.get("link", "")
            if not fa_link:
                fa_link = "FA"
            fa_seg_str = (
                f'{fa_link}{fa_matched_dct.get("c", 0)}:{fa_matched_dct.get("db", 0)}'
            )
            if fa_seg_str.lower().startswith("none"):
                fa_seg_str = fa_seg_str[4:]
            if mod_code:
                mod_obj = Modifications(mod_code)
                fa_info_dct["mod_obj"] = mod_obj
                fa_info_dct["mod_info"] = mod_obj.sum_mod_info
                fa_info_dct["mod_id"] = mod_obj.mod_id
                fa_info_dct["fa_level"] = f"{mod_obj.mod_level}"
                fa_linked_ids_dct = {}
                mod_linked_ids_dct = mod_obj.mod_linked_ids  # type: dict
                for lv in mod_linked_ids_dct:
                    fa_linked_ids_dct[f"{lv}"] = f"{fa_seg_str}{mod_linked_ids_dct[lv]}"
                fa_info_dct["fa_linked_ids"] = fa_linked_ids_dct

                if mod_obj.mod_list:
                    for mod_seg in mod_obj.mod_list:
                        if (
                            mod_seg.get("cv", "DB") != "DB"
                            and mod_seg.get("count", 0) != 0
                        ):
                            is_modified = True
            else:
                if fa_matched_dct.get("db", 0) == 0:
                    fa_info_dct["fa_level"] = "4.2"
                    fa_info_dct["fa_linked_ids"] = {
                        "0": fa_seg_str,
                        "1": fa_seg_str,
                        "2": fa_seg_str,
                        "3": fa_seg_str,
                        "3.1": fa_seg_str,
                        "3.2": fa_seg_str,
                        "4": fa_seg_str,
                        "4.1": fa_seg_str,
                        "4.2": fa_seg_str,
                    }
                else:
                    fa_info_dct["fa_level"] = "4"
                    fa_info_dct["fa_linked_ids"] = {
                        "0": fa_seg_str,
                        "1": fa_seg_str,
                        "2": fa_seg_str,
                        "3": fa_seg_str,
                        "4": fa_seg_str,
                    }
                if self.lipid_code != fa_seg_str:
                    self.lipid_code = fa_seg_str
            fa_info_dct["fa_info"] = {
                "c": int(fa_matched_dct.get("c", 0)),
                "db": int(fa_matched_dct.get("db", 0)),
                "link": fa_link,
                "is_modified": is_modified,
            }
            return fa_info_dct
        else:
            raise ValueError(f"Cannot parse FA sting: {self.lipid_code}")

    def to_json(self):
        fa_lite_info_dct = self.fa_info_dct
        fa_lite_info_dct.pop("mod_obj", None)
        fa_json_str = json.dumps(fa_lite_info_dct)
        if check_json(validator=self.validator, json_obj=json.loads(fa_json_str)):
            return fa_json_str
        else:
            raise Exception(f"JSON Schema check FAILED. Schema {self.schema}")


if __name__ == "__main__":

    hg_lst = [
        "PA",
        "PC",
        "PE",
        "PG",
        "PI",
        "PS",
        "SM",
        "SPB",
        "Cer",
        "PIP",
        "PIP2",
        "PIP3",
    ]

    for hg in hg_lst:
        hg_obj = HeadGroup(hg_code=hg)
        logger.info(hg_obj.to_json())

    mod_code_lst = [
        r"16:0",
        r"18:1",
        r"18:1<{9Z}>",
        r"20:4<-18>",
        r"20:4<+46>",
        r"20:4<+3O,-2H>",
        r"20:4<2OH,Ke>",
        r"20:4<2OH{8,11},Ke{14}>",
        r"20:4<{5,9,12,15},2OH{8,11},Ke{14}>",
        r"20:4<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>",
        r"20:4<2OH{8R,11S},Ke{14}>",
        r"20:4<{5,9,12,15},2OH{8R,11S},Ke{14}>",
        r"20:4<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>",
        r"FA16:0",
        r"FA18:1",
        r"FA18:1<{9Z}>",
        r"FA20:4<-18>",
        r"FA20:4<+46>",
        r"FA20:4<+3O,-2H>",
        r"FA20:4<2OH,Ke>",
        r"FA20:4<2OH{8,11},Ke{14}>",
        r"FA20:4<{5,9,12,15},2OH{8,11},Ke{14}>",
        r"FA20:4<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>",
        r"FA20:4<2OH{8R,11S},Ke{14}>",
        r"FA20:4<{5,9,12,15},2OH{8R,11S},Ke{14}>",
        r"FA20:4<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>",
        r"O-16:0",
        r"O-18:1",
        r"O-18:1<{9Z}>",
        r"O-20:4<-18>",
        r"O-20:4<+46>",
        r"O-20:4<+3O,-2H>",
        r"O-20:4<2OH,Ke>",
        r"O-20:4<2OH{8,11},Ke{14}>",
        r"O-20:4<{5,9,12,15},2OH{8,11},Ke{14}>",
        r"O-20:4<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>",
        r"O-20:4<2OH{8R,11S},Ke{14}>",
        r"O-20:4<{5,9,12,15},2OH{8R,11S},Ke{14}>",
        r"O-20:4<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>",
        r"P-16:0",
        r"P-18:1",
        r"P-18:1<{9Z}>",
        r"P-20:4<-18>",
        r"P-20:4<+46>",
        r"P-20:4<+3O,-2H>",
        r"P-20:4<2OH,Ke>",
        r"P-20:4<2OH{8,11},Ke{14}>",
        r"P-20:4<{5,9,12,15},2OH{8,11},Ke{14}>",
        r"P-20:4<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>",
        r"P-20:4<2OH{8R,11S},Ke{14}>",
        r"P-20:4<{5,9,12,15},2OH{8R,11S},Ke{14}>",
        r"P-20:4<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>",
    ]

    for usr_mod_code in mod_code_lst:
        logger.info(f"Test FA str: {usr_mod_code}")

        fa_obj = FattyAcid(usr_mod_code)

        # logger.debug(f"to all levels: {fa_obj.to_all_levels()}")
        # logger.debug(
        #     f"to all levels without <> : {fa_obj.to_all_levels(angle_brackets=False)}"
        # )
        #
        fa_json = fa_obj.to_json()
        logger.debug(fa_json)
        fa_ids_dct = json.loads(fa_json).get("fa_linked_ids")
        logger.info("".join([f"\n {s}: {fa_ids_dct[s]}" for s in fa_ids_dct]))

    logger.info("FINISHED")
