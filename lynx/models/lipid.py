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

import itertools
import json
from operator import itemgetter
import re

from jsonschema import Draft7Validator
from natsort import natsorted

from ..models.defaults import (
    lynx_schema_cfg,
    lipid_level_lst,
    mod_db_level_lst,
)
from lynx.utils.basics import get_abs_path
from lynx.utils.cfg_reader import api_version
from lynx.utils.log import app_logger
from lynx.utils.toolbox import check_json


class Lipid(object):
    def __init__(self, lipid_code: str, logger=app_logger):

        self.lipid_code = lipid_code
        self.lynx_class_lv0 = ""
        self.schema = "lynx_core"
        with open(get_abs_path(lynx_schema_cfg[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(json.load(s_obj))

        self.level = "B0"
        self._lipid_level = "B"
        self._max_mod_level = 0
        self.is_modified = False
        self.sum_info = self.__post_init__()

        self.residues = self.sum_info.get("residues", [])
        self.level = self.sum_info.get("level", "")
        self.linked_ids = self.sum_info.get("linked_ids", {})
        self.logger = logger
        self.logger.info(
            f"Level {self.level:4s} FattyAcid created from: {self.lipid_code}"
        )

    def __identify_class__(self):
        lipid_class = ""
        lipid_segments = {}
        fa_match = fa_rgx.match(self.lipid_code)
        pl_match = pl_rgx.match(self.lipid_code)
        sp_match = sp_rgx.match(self.lipid_code)
        gl_match = gl_rgx.match(self.lipid_code)
        cl_match = cl_rgx.match(self.lipid_code)

        if fa_match:
            lipid_class = "FA"
            lipid_segments = fa_match.groupdict()
        if pl_match:
            lipid_segments = pl_match.groupdict()
            lipid_class = lipid_segments.get("hg_class", "")
        if sp_match:
            lipid_segments = sp_match.groupdict()
            lipid_class = lipid_segments.get("hg_class", "")
        if gl_match:
            lipid_segments = gl_match.groupdict()
            fa_list = [
                lipid_segments.get("fa1", None),
                lipid_segments.get("fa2", None),
                lipid_segments.get("fa3", None),
            ]
            fa_list = [f for f in fa_list if f and len(f) > 2]
            if len(fa_list) == 3:
                lipid_class = "TG"
            elif len(fa_list) == 2:
                lipid_class = "DG"
            elif len(fa_list) == 1:
                lipid_class = "MG"
        if cl_match:
            lipid_class = "CL"
            lipid_segments = fa_match.groupdict()
        if lipid_class:
            self.lynx_class_lv0 = lipid_class
        return lipid_segments

    def __identify_fa__(self, lipid_segments, fa_count: int = 1):
        fa_info_dct: dict = {}
        fa = f"fa{fa_count}"
        sn = f"sn{fa_count}"
        if fa_count > 1:
            ps = f"position{fa_count-1}"
        else:
            ps = f"position1"
        if fa_count > 0:
            if fa in lipid_segments and lipid_segments.get(fa, None):
                if lipid_segments.get(ps, None) == "/":
                    self._lipid_level = "S"
                    fa_info_dct = {"position": sn, "id": lipid_segments[fa]}
                else:
                    if lipid_segments.get(ps, None) == "_":
                        self._lipid_level = "D"
                    fa_info_dct = {"id": lipid_segments[fa]}
        return fa_info_dct

    def __identify_level__(self, lipid_segments: dict):

        max_mod_level = 0
        max_unmod_level = 0
        if self.lynx_class_lv0 == "FA" and "link" in lipid_segments:
            res_info_lst = [{"id": self.lipid_code}]
            self._lipid_level = "S"
        else:
            res_info_lst = [
                self.__identify_fa__(lipid_segments, i)
                for i in [1, 2, 3, 4, 5]
                if self.__identify_fa__(lipid_segments, i)
            ]
        if self._lipid_level != "S":  # Sort all residues by names in B and D level
            res_info_lst = natsorted(res_info_lst, key=itemgetter(*["id"]))

        if (
            "hg_class" in lipid_segments
            and lipid_segments.get("hg_class", None)
            and "fa3" not in lipid_segments
        ):
            if lipid_segments.get("position1", None) == "/":
                res_info_lst.append(
                    {"position": "sn3", "id": lipid_segments["hg_class"]}
                )
            else:
                res_info_lst.append({"id": lipid_segments["hg_class"]})

        for res in res_info_lst:
            res_type = ""
            if fa_rgx.match(res["id"]):
                res_type = "FA"
                res_obj = FattyAcid(res["id"])
                res_level = res_obj.fa_level
                if res_obj.is_modified:
                    max_mod_level = max(max_mod_level, float(res_level))
                    self.is_modified = True
                else:
                    max_unmod_level = max(max_unmod_level, float(res_level))
            else:
                try:
                    res_obj = LipidClass(res["id"])
                    res_type = "HG"
                except Exception as err:
                    self.logger.error(err)
                    raise ValueError(f'Cannot parse string as HeadGroup: {res["id"]}')
            res["info"] = json.loads(res_obj.to_json())
            res["info"]["type"] = res_type

        if max_unmod_level > max_mod_level:
            max_unmod_level = max_mod_level
        self._max_mod_level = max_mod_level
        self.level = f"{self._lipid_level}{max_mod_level}"
        if self.level.endswith(".0"):
            self.level = self.level[:-2]

        return res_info_lst

    def _get_linked_ids(self, res_lst):

        linked_ids_dct = {}

        other_pre_dct = {}
        lv_lst = self.get_levels()
        fa_res_dct = {}
        fa_info_lst = []
        for mod_lv in lv_lst["mod_lv_lst"]:
            fa_res_dct[mod_lv] = []
        self.logger.info(res_lst)
        for res in res_lst:
            res_info = res["info"]
            if res_info.get("type", None) == "FA":
                res_obj = FattyAcid(res_info.get("id", None))
                res_level = res_obj.fa_level
                res_id = res_obj.id
                lift_res_level = 0
                fa_linked_ids = res_obj.fa_linked_ids
                if not res_obj.is_modified and float(res_level) < self._max_mod_level:
                    lift_res_level = self._max_mod_level - float(res_level)
                res_id_dct = {}
                if self._max_mod_level < float(res_level):
                    for lv in fa_linked_ids:
                        if float(lv) <= self._max_mod_level:
                            res_id_dct[lv] = fa_linked_ids[lv]
                else:
                    res_id_dct = fa_linked_ids
                if res_obj.is_modified:
                    for mod_lv in fa_res_dct:
                        fa_res_dct[mod_lv].append(res_id_dct[mod_lv].strip("FA"))
                else:
                    for mod_lv in fa_res_dct:
                        if len(mod_lv) == 1:
                            fa_res_dct[mod_lv].append(
                                res_id_dct.get("0", "").strip("FA")
                            )
                        elif mod_lv.endswith(".1"):
                            if res_id_dct.get("0.1", "").strip("FA"):
                                fa_res_dct[mod_lv].append(
                                    res_id_dct.get("0.1", "").strip("FA")
                                )
                            else:
                                fa_res_dct[mod_lv].append(
                                    res_id_dct.get("0", "").strip("FA")
                                )

                        elif mod_lv.endswith(".2"):
                            if res_id_dct.get("0.2", "").strip("FA"):
                                fa_res_dct[mod_lv].append(
                                    res_id_dct.get("0.2", "").strip("FA")
                                )
                            elif res_id_dct.get("0.1", "").strip("FA"):
                                fa_res_dct[mod_lv].append(
                                    res_id_dct.get("0.1", "").strip("FA")
                                )
                            else:
                                fa_res_dct[mod_lv].append(
                                    res_id_dct.get("0", "").strip("FA")
                                )

                if float(res_level) > 3:
                    bulk_level = "3"
                else:
                    bulk_level = res_level
                fa_info_lst.append(res_obj.to_segments(mod_level=bulk_level))
            elif res_info.get("type", None) == "HG":
                res_obj = LipidClass(res_info.get("id", None))
                other_pre_dct["HG"] = res_obj.id
        bulk_linked_ids = {}
        if lv_lst["lynx_lv_lst"][0].startswith("B"):
            sum_c = 0
            sum_db = 0
            sum_link = []
            sum_mod_lst = []
            for info in fa_info_lst:
                sum_c += int(info["c"])
                sum_db += int(info["db"])
                sum_link.append(info["link"].strip("FA"))
                sum_mod_lst.append(info["mod_id"].strip("<>"))
            sum_fa_str = f"{''.join(sum_link)}{sum_c}:{sum_db}<{','.join(sum_mod_lst).strip(',')}>"
            bulk_linked_ids = FattyAcid(sum_fa_str).fa_linked_ids

        for lv in lv_lst["lynx_lv_lst"]:
            if self.lynx_class_lv0 == "FA":
                fa_str_lst = fa_res_dct.get(lv[1:], [])
                if len(fa_str_lst) == 1:
                    fa_str = fa_str_lst[0]
                    if fa_str[:2] in ["FA", "O-", "P-"]:
                        pass
                    else:
                        fa_str = f"FA{fa_str}"
                    linked_ids_dct[lv] = fa_str
            else:
                if lv[0] == "S":
                    fa_seg_str = "/".join(fa_res_dct.get(lv[1:], []))
                    if "HG" in other_pre_dct:
                        linked_ids_dct[lv] = f'{other_pre_dct["HG"]}({fa_seg_str})'
                    else:
                        linked_ids_dct[lv] = f"{fa_seg_str}"
                elif lv[0] == "D":
                    pre_fa_seg_lst = fa_res_dct.get(lv[1:], [])
                    fa_seg_lst = []
                    # sort O-/P- first in D level
                    for fa_seg in pre_fa_seg_lst:
                        if fa_seg.startswith("O-"):
                            fa_seg_lst.append(fa_seg)
                            pre_fa_seg_lst.remove(fa_seg)
                        elif fa_seg.startswith("P-"):
                            fa_seg_lst.append(fa_seg)
                            pre_fa_seg_lst.remove(fa_seg)
                    fa_seg_lst = natsorted(fa_seg_lst) + natsorted(pre_fa_seg_lst)
                    fa_seg_str = "_".join(fa_seg_lst)
                    if "HG" in other_pre_dct:
                        linked_ids_dct[lv] = f'{other_pre_dct["HG"]}({fa_seg_str})'
                    else:
                        linked_ids_dct[lv] = f"{fa_seg_str}"
                elif lv[0] == "B" and int(lv[1]) <= 3:
                    fa_seg_str = bulk_linked_ids.get(lv[1:], "").strip("FA")
                    if "HG" in other_pre_dct:
                        linked_ids_dct[lv] = f'{other_pre_dct["HG"]}({fa_seg_str})'
                    else:
                        linked_ids_dct[lv] = f"{fa_seg_str}"
        return linked_ids_dct

    def __post_init__(self):
        lipid_segments = self.__identify_class__()

        if self.lynx_class_lv0 and lipid_segments:
            res_info_lst = self.__identify_level__(lipid_segments)
            lipid_info_dct = {
                "api_version": api_version,
                "id": self.lipid_code,
                "level": self.level,
                "linked_ids": self._get_linked_ids(res_info_lst),
                "is_modified": self.is_modified,
                "info": {"main_class": self.lynx_class_lv0, "residues": res_info_lst},
            }
            self.logger.info(f"modification level: {self._max_mod_level}")
            self.logger.info(f"\n{lipid_info_dct}")

            return lipid_info_dct
        else:
            raise ValueError(f"Cannot parse Lipid sting: {self.lipid_code}")

    def to_json(self):
        sum_info = self.sum_info
        sum_info.pop("mod_obj", None)
        json_str = json.dumps(sum_info)
        if check_json(
            validator=self.validator, json_obj=json.loads(json_str, logger=self.logger)
        ):
            return json_str
        else:
            raise Exception(f"JSON Schema check FAILED. Schema {self.schema}")

    def get_levels(self):
        out_lipid_levels_lst = []
        if self._lipid_level in lipid_level_lst:
            l_lv_idx = lipid_level_lst.index(self._lipid_level)
            out_lipid_levels_lst = lipid_level_lst[: l_lv_idx + 1]

        max_level_str = str(self._max_mod_level)
        if max_level_str.endswith(".0"):
            max_level_str = max_level_str[:-2]
        if max_level_str in mod_db_level_lst:
            mod_idx = mod_db_level_lst.index(max_level_str)
            out_mod_levels_lst = mod_db_level_lst[: mod_idx + 1]
        else:
            raise ValueError(f"Modification level not supported: {max_level_str}")
        if len(max_level_str) == 1:
            out_mod_levels_lst = [
                out_lv for out_lv in out_mod_levels_lst if len(out_lv) == 1
            ]
        elif max_level_str.endswith(".1"):
            out_mod_levels_lst = [
                out_lv
                for out_lv in out_mod_levels_lst
                if len(out_lv) == 1 or out_lv.endswith(".1")
            ]
        pre_lynx_lv_lst = [
            "".join(s)
            for s in itertools.product(out_lipid_levels_lst, out_mod_levels_lst)
        ]
        lynx_lv_lst = []
        for lv in pre_lynx_lv_lst:
            if not re.match(r"B\d\.\d", lv):
                lynx_lv_lst.append(lv)
        lv_info_dct = {
            "lipid_lv_lst": out_lipid_levels_lst,
            "mod_lv_lst": out_mod_levels_lst,
            "lynx_lv_lst": lynx_lv_lst,
        }

        return lv_info_dct


if __name__ == "__main__":

    hg_lst = [
        "PA",
        "PC",
        "PE",
        "PG",
        "PI",
        "PS",
        "SM",
        # "SPB",
        # "Cer",
        "PIP",
        "PIP2",
        "PIP3",
    ]

    sn1_lst = ["18:1<{9Z}>", "18:1", "16:0", "O-16:0", "P-16:0", "O-18:1", "P-18:1"]
    # sn1_lst = ["18:1<{9Z}>"]

    sn2_lst = [
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
        # r"O-20:4<-18>",
        # r"O-20:4<+46>",
        # r"O-20:4<+3O,-2H>",
        # r"O-20:4<2OH,Ke>",
        # r"O-20:4<2OH{8,11},Ke{14}>",
        # r"O-20:4<{5,9,12,15},2OH{8,11},Ke{14}>",
        # r"O-20:4<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>",
        # r"O-20:4<2OH{8R,11S},Ke{14}>",
        # r"O-20:4<{5,9,12,15},2OH{8R,11S},Ke{14}>",
        # r"O-20:4<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>",
        # r"P-20:4<-18>",
        # r"P-20:4<+46>",
        # r"P-20:4<+3O,-2H>",
        # r"P-20:4<2OH,Ke>",
        # r"P-20:4<2OH{8,11},Ke{14}>",
        # r"P-20:4<{5,9,12,15},2OH{8,11},Ke{14}>",
        # r"P-20:4<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>",
        # r"P-20:4<2OH{8R,11S},Ke{14}>",
        # r"P-20:4<{5,9,12,15},2OH{8R,11S},Ke{14}>",
        # r"P-20:4<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>",
    ]

    pl_d_lst = []
    pl_s_lst = []

    for hg in hg_lst:
        for sn1 in sn1_lst:
            for sn2 in sn2_lst:
                if "-" not in sn1 and "-" not in sn2:
                    pl_d_lst.append(f"{hg}({sn1}_{sn2})")
                    pl_s_lst.append(f"{hg}({sn1}/{sn2})")

    # lipid_lst = sn1_lst + sn2_lst + pl_s_lst + pl_d_lst
    lipid_lst = sn1_lst
    counter = 0
    for pl in lipid_lst:
        counter += 1
        self.logger.info(f"Test Lipid #{counter} : {pl}")
        pl_obj = Lipid(lipid_code=pl)
        self.logger.info(f"Export JSON \n {pl_obj.to_json()}")
        self.logger.info(
            "".join([f"\n {s}: {pl_obj.linked_ids[s]}" for s in pl_obj.linked_ids])
        )

    self.logger.info("FINISHED")
