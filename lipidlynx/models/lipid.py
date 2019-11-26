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
from lipidlynx.models.patterns import fa_rgx, pl_rgx, sp_rgx, gl_rgx, cl_rgx
from lipidlynx.models.residues import HeadGroup, FattyAcid


class Lipid(object):
    def __init__(self, lipid_code: str):

        self.lipid_code = lipid_code
        self.schema = "lynx_core"
        with open(get_abs_path(lynx_schema[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(json.load(s_obj))

        self.sum_info = self.__post_init__()

        self.residues = self.sum_info.get("residues", [])
        self.level = self.sum_info.get("level", "")
        self.linked_ids = self.sum_info.get("linked_ids", {})
        logger.info(f"Level {self.level:4s} FattyAcid created from: {self.lipid_code}")

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

        return lipid_class, lipid_segments

    def __post_init__(self):
        ident_class, lipid_segments = self.__identify_class__()
        res_info_lst = []
        res_obj_lst = []
        lipid_level = 'B'
        max_mod_level = 0
        max_unmod_level = 0
        linked_ids_dct = {}

        if ident_class and lipid_segments:
            if "fa1" in lipid_segments and lipid_segments.get("fa1", None):
                if lipid_segments.get("position1", None) == '/':
                    lipid_level = 'S'
                    res_info_lst.append(
                        {"position": "sn1", "residue_id": lipid_segments["fa1"]}
                    )
                else:
                    if lipid_segments.get("position1", None) == '_':
                        lipid_level = 'D'
                    res_info_lst.append({"residue_id": lipid_segments["fa1"]})
                res_obj_lst.append(FattyAcid(lipid_segments["fa1"]))
            if "fa2" in lipid_segments and lipid_segments.get("fa2", None):
                if lipid_segments.get("position1", None) == '/':
                    lipid_level = 'S'
                    res_info_lst.append(
                        {"position": "sn2", "residue_id": lipid_segments["fa2"]}
                    )
                else:
                    if lipid_segments.get("position1", None) == '_':
                        lipid_level = 'D'
                    res_info_lst.append({"residue_id": lipid_segments["fa2"]})
                res_obj_lst.append(FattyAcid(lipid_segments["fa2"]))
            if "fa3" in lipid_segments and lipid_segments.get("fa3", None):
                if lipid_segments.get("position2", None) == '/':
                    lipid_level = 'S'
                    res_info_lst.append(
                        {"position": "sn3", "residue_id": lipid_segments["fa3"]}
                    )
                else:
                    if lipid_segments.get("position2", None) == '_':
                        lipid_level = 'D'
                    res_info_lst.append({"residue_id": lipid_segments["fa3"]})
                res_obj_lst.append(FattyAcid(lipid_segments["fa3"]))
            if "fa4" in lipid_segments and lipid_segments.get("fa4", None):
                if lipid_segments.get("position3", None) == '/':
                    lipid_level = 'S'
                    res_info_lst.append(
                        {"position": "sn4", "residue_id": lipid_segments["fa4"]}
                    )
                else:
                    if lipid_segments.get("position3", None) == '_':
                        lipid_level = 'D'
                    res_info_lst.append({"residue_id": lipid_segments["fa4"]})
                res_obj_lst.append(FattyAcid(lipid_segments["fa4"]))
            if (
                "hg_class" in lipid_segments
                and lipid_segments.get("hg_class", None)
                and "fa3" not in lipid_segments
            ):
                if lipid_segments.get("position1", None) == '/':
                    res_info_lst.append(
                        {"position": "sn3", "residue_id": lipid_segments["hg_class"]}
                    )
                else:
                    res_info_lst.append({"residue_id": lipid_segments["hg_class"]})

            for res_obj in res_obj_lst:
                if isinstance(res_obj, FattyAcid):
                    res_level = res_obj.fa_level
                    if res_obj.is_modified:
                        max_mod_level = max(max_mod_level, float(res_level))
                    else:
                        max_unmod_level = max(max_unmod_level, float(res_level))

            if max_unmod_level > max_mod_level:
                max_unmod_level = max_mod_level
            lynx_level = f'{lipid_level}{max_mod_level}'
            if lynx_level.endswith('.0'):
                lynx_level = lynx_level[:-2]
            linked_ids_dct[lynx_level] = self.lipid_code
            lipid_info_dct = {
                "_version": lynx_schema.get("_version", "0.1"),
                "id": self.lipid_code,
                "lipid_class": {'main_class': ident_class},
                "level": lynx_level,
                "residues": res_info_lst,
                "linked_ids": linked_ids_dct
            }

            logger.info(f'modification level: {max_mod_level}')
            logger.info(f'\n{lipid_info_dct}')

            return lipid_info_dct
        else:
            raise ValueError(f"Cannot parse Lipid sting: {self.lipid_code}")

    def to_json(self):
        sum_info = self.sum_info
        sum_info.pop("mod_obj", None)
        json_str = json.dumps(sum_info)
        if check_json(validator=self.validator, json_obj=json.loads(json_str)):
            return json_str
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

    sn1_lst = ["16:0", "O-16:0", "P-16:0", "18:1", "O-18:1", "P-18:1"]

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

    pl_d_lst = []
    pl_s_lst = []

    for hg in hg_lst:
        for sn1 in sn1_lst:
            for sn2 in sn2_lst:
                if '-' not in sn1 and '-' not in sn2:
                    pl_d_lst.append(f"{hg}({sn1}_{sn2})")
                    pl_s_lst.append(f"{hg}({sn1}/{sn2})")

    pl_lst = pl_d_lst + pl_s_lst
    counter = 0
    for pl in pl_lst:
        counter += 1
        logger.info(f"Test Lipid #{counter} : {pl}")
        pl_obj = Lipid(lipid_code=pl)
        logger.info(pl_obj.to_json())

    logger.info("FINISHED")
