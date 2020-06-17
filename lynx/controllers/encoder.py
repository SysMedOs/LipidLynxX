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

from typing import Dict, List

from natsort import natsorted

from lynx.utils.params_loader import load_output_rule
from lynx.controllers.decoder import Decoder
from lynx.models.residue import Residue, merge_residues
from lynx.models.defaults import (
    default_output_rules,
    default_input_rules,
    supported_levels,
)
from lynx.utils.log import app_logger


class Encoder(object):
    def __init__(
        self,
        output_rules: dict = default_output_rules,
        style: str = "LipidLynxX",
        input_rules: dict = default_input_rules,
        logger=app_logger,
    ):
        self.export_rule = style
        self.output_rules = load_output_rule(output_rules, style)
        self.class_rules = self.output_rules.get("LMSD_CLASSES", {})
        self.mod_rules = self.output_rules.get("MODS", {}).get("MOD", {})
        self.sum_mod_rules = self.output_rules.get("MODS", {}).get("SUM_MODS", {})
        self.residue_rules = self.output_rules.get("RESIDUES", {}).get("RESIDUE", {})
        self.separator_levels = self.output_rules.get("SEPARATORS", {}).get(
            "SEPARATOR_LEVELS", {}
        )
        self.separators = self.output_rules.get("SEPARATORS", {})
        self.extractor = Decoder(rules=input_rules)
        self.logger = logger

    @staticmethod
    def get_best_id(candidate: Dict[str, str]) -> str:
        c_lv_lst = natsorted(list(candidate.keys()))
        c_max_str = candidate.get(c_lv_lst[-1], "")
        return c_max_str

    @staticmethod
    def get_best_id_series(candidates: List[dict]) -> dict:
        best_id_dct = {}
        num_lv = 0
        max_str = ""
        for c_info in candidates:
            for c in c_info:
                c_lv_lst = natsorted(list(c_info[c].keys()))
                c_num_lv = len(c_lv_lst)
                c_max_str = c_info[c].get(c_lv_lst[-1], "")
                if c_num_lv > num_lv:
                    best_id_dct = c_info[c]
                    max_str = c_max_str
                    num_lv = c_num_lv
                else:
                    c_max_str = c_info[c].get(c_lv_lst[-1], "")
                    if len(c_max_str) > len(max_str):
                        best_id_dct = c_info[c]
                        max_str = c_max_str
                        num_lv = c_num_lv
                    else:
                        pass

        # add levels for B0, D0, S0 lipids
        if best_id_dct:
            if list(best_id_dct.keys())[-1].endswith("0"):
                for lv in list(best_id_dct.keys()):
                    for lv_base in ["B", "D", "S"]:
                        if lv.startswith(lv_base):
                            lv_id = best_id_dct.get(lv)
                            for add_lv in supported_levels:
                                if add_lv.startswith(lv_base) and add_lv != lv:
                                    best_id_dct[add_lv] = lv_id

        if best_id_dct.get("D1") and not best_id_dct.get("B1"):
            best_id_dct["B1"] = best_id_dct["D1"]
        if best_id_dct.get("D2") and not best_id_dct.get("B2"):
            best_id_dct["B2"] = best_id_dct["D2"]

        return best_id_dct

    # def check_rest(self, segment_text: str, segment_name: str, lmsd_class: str):
    #     patterns_dct = self.class_rules[lmsd_class].get(segment_name)
    #     out_seg_lst = []
    #     if segment_text and patterns_dct:
    #         for s_rgx in patterns_dct:
    #             self.logger.debug(
    #                 f"Test {segment_text} on {segment_name} of {lmsd_class} using {s_rgx}"
    #             )
    #             s_matched = s_rgx.match(segment_text)
    #             if s_matched:
    #                 defined_seg = patterns_dct[s_rgx]
    #                 if defined_seg == "EXCEPTIONS":
    #                     if lmsd_class in ["GP12"]:
    #                         out_seg_lst.append(segment_text)
    #                     if lmsd_class in ["SP05", "SP06"]:
    #                         out_seg_lst.append(segment_text)
    #                 else:
    #                     out_seg_lst.append(defined_seg)
    #
    #         out_seg_lst = list(filter(None, list(set(out_seg_lst))))
    #
    #     return self.get_best_candidate(out_seg_lst)

    def get_residues(self, residues: dict):
        residues_order = residues.get("RESIDUES_ORDER", [])
        residues_sep_level = residues.get("RESIDUES_SEPARATOR_LEVEL", "")
        residues_info = residues.get("RESIDUES_INFO", [])
        res_count = len(residues_order)
        sum_residues_str = ""
        res_lv_id_dct = {}
        res_lv_dct = {}
        sum_lv_lst = []
        for res_abbr in residues_info:
            res_obj = Residue(residues_info[res_abbr], nomenclature=self.export_rule)
            res_lv_id_dct[res_abbr] = res_obj.linked_ids
            res_lv_dct[res_abbr] = list(res_obj.linked_ids.keys())
            sum_lv_lst.extend(res_lv_dct[res_abbr])

        sum_lv_lst = natsorted(set(sum_lv_lst))
        for res in res_lv_id_dct:
            r_lv_dct = res_lv_id_dct[res]
            if list(r_lv_dct.keys()) == ["0"]:
                for lv in sum_lv_lst:
                    r_lv_dct[lv] = r_lv_dct["0"]
        res_lv_id_lst_dct = {}
        for sum_lv in sum_lv_lst:
            lv_id_lst = []
            for r in residues_order:
                r_lv_id = res_lv_id_dct.get(r, {}).get(sum_lv, None)
                if r_lv_id:
                    lv_id_lst.append(r_lv_id)
            if len(lv_id_lst) == res_count:
                res_lv_id_lst_dct[sum_lv] = lv_id_lst

        sum_res_id_lv_dct = {}
        sum_res_sep_lv_lst = []

        # set FA or class with one Res into level s
        if len(residues_order) == 1 and len(sum_lv_lst) > 0:
            if "0.1" in sum_lv_lst and residues_sep_level == "B":
                residues_sep_level = "S"

        if residues_sep_level == "S":
            sum_res_sep_lv_lst = ["B", "D", "S"]
        elif residues_sep_level == "D":
            sum_res_sep_lv_lst = ["B", "D"]
        elif residues_sep_level == "B":
            sum_res_sep_lv_lst = ["B"]

        for sep_lv in sum_res_sep_lv_lst:
            if sep_lv == "B":
                # prepare bulk level
                merged_res_obj = merge_residues(
                    residues_info, nomenclature=self.export_rule
                )
                merged_res_linked_ids = merged_res_obj.linked_ids
                merged_res_lv_lst = list(merged_res_obj.linked_ids.keys())
                for merged_res_lv in merged_res_lv_lst:
                    sum_res_id_lv_dct[f"B{merged_res_lv}"] = merged_res_linked_ids[
                        merged_res_lv
                    ]
            else:
                for res_lv in res_lv_id_lst_dct:
                    sum_res_id_lv_dct[
                        f"{sep_lv}{res_lv}"
                    ] = f"{self.separator_levels[sep_lv]}".join(
                        res_lv_id_lst_dct[res_lv]
                    )

        return sum_res_id_lv_dct

    def check_segments(self, parsed_info: dict):
        segments_dct = {}
        lmsd_classes = parsed_info.get("LMSD_CLASSES", None)
        segments = parsed_info["SEGMENTS"]
        c_prefix_lst = segments.get("PREFIX", [])
        c_suffix_lst = segments.get("SUFFIX", [])
        c_has_prefix = False
        if len(c_prefix_lst) == 1 and c_prefix_lst[0]:
            c_prefix_seg = c_prefix_lst[0]
            c_has_prefix = True
        else:
            c_prefix_seg = ""
        c_has_suffix = False
        if len(c_suffix_lst) == 1 and c_suffix_lst[0]:
            c_suffix_seg = c_suffix_lst[0]
            c_has_suffix = True
        else:
            c_suffix_seg = ""
        residues = parsed_info.get("RESIDUES", {})
        sum_res_id_lv_dct = self.get_residues(residues)
        obs_c_seg_lst = segments.get("CLASS", [])
        c_seg = ""
        if obs_c_seg_lst and len(obs_c_seg_lst) == 1:
            obs_c_seg = obs_c_seg_lst[0]
            for c in lmsd_classes:
                c_segments_dct = {}
                c_orders = []
                if c in self.class_rules:
                    c_orders = self.class_rules[c].get("ORDER", [])
                    c_identifier_dct = self.class_rules[c].get("CLASS", {})
                    for c_identifier in c_identifier_dct:
                        if c_identifier.match(obs_c_seg):
                            c_seg = c_identifier_dct.get(c_identifier, "")
                            for lv in sum_res_id_lv_dct:
                                c_lv_segments = {
                                    "CLASS": c_seg,
                                    "SUM_RESIDUES": sum_res_id_lv_dct[lv],
                                }
                                if c_has_prefix:
                                    c_lv_segments["PREFIX"] = c_prefix_seg
                                if c_has_suffix:
                                    c_lv_segments["SUFFIX"] = c_suffix_seg
                                c_segments_dct[lv] = c_lv_segments
                        else:
                            pass
                else:
                    pass
                if c_segments_dct and c_orders:
                    segments_dct[c] = {"ORDER": c_orders, "INFO": c_segments_dct}
        else:
            self.logger.warning(f"No Class identified!")

        return segments_dct

    def compile_segments(self, segments: dict):
        comp_seg_dct = {}
        for c in segments:
            c_lv_id_dct = {}
            c_seg_dct = segments[c]
            c_seg_order = c_seg_dct.get("ORDER", [])
            c_comp_seg_dct = c_seg_dct.get("INFO", {})
            c_optional_seg = self.class_rules.get(c, {}).get("OPTIONAL", [])
            for lv in c_comp_seg_dct:
                lv_seg_info = c_comp_seg_dct[lv]
                lv_seg_lst = []
                for c_seg in c_seg_order:
                    if c_seg in lv_seg_info:
                        lv_seg_lst.append(lv_seg_info[c_seg])
                    elif c_seg in self.separators and c_seg != "SEPARATOR_LEVELS":
                        lv_seg_lst.append(self.separators[c_seg])
                    else:
                        if c_seg not in c_optional_seg:
                            self.logger.debug(
                                f"Segments not found: {c_seg}, defined orders {c_seg_order}"
                            )
                        else:
                            pass
                c_lv_id_dct[lv] = "".join(lv_seg_lst)
            comp_seg_dct[c] = c_lv_id_dct

        return comp_seg_dct

    def export_all_levels(self, lipid_name: str) -> dict:

        extracted_info = self.extractor.extract(lipid_name)
        export_info = []
        if extracted_info:
            for p in extracted_info:
                p_info = extracted_info[p]
                self.logger.info(p_info)
                for in_r in p_info:
                    r_info = p_info[in_r]  # type: dict
                    checked_seg_info = self.check_segments(r_info)
                    comp_dct = self.compile_segments(checked_seg_info)
                    export_info.append(comp_dct)
            best_export_dct = self.get_best_id_series(export_info)
            self.logger.debug(f"Convert Lipid: {lipid_name} into:\n{best_export_dct}")
        else:
            best_export_dct = {}

        return best_export_dct

    def convert(self, lipid_name: str, level: str = None) -> str:
        if level in supported_levels:
            best_id = self.export_level(lipid_name, level=level)
        else:
            all_lv_id_dct = self.export_all_levels(lipid_name)
            best_id = ""
            if all_lv_id_dct:
                best_id = self.get_best_id(all_lv_id_dct)
            else:
                pass
        return best_id

    def export_level(
        self, lipid_name: str, level: str = "B1",
    ):

        lv_id = ""
        all_lv_id_dct = self.export_all_levels(lipid_name)
        if level in supported_levels:
            if level in all_lv_id_dct:
                lv_id = all_lv_id_dct[level]
            else:
                self.logger.warning(
                    f"Lipid: {lipid_name} cannot be converted into level: {level}. "
                    f"Can be converted into: {all_lv_id_dct}"
                )
        else:
            self.logger.warning(
                f"Level: {level} not supported. Supported levels: {supported_levels}"
            )

        return lv_id

    def export_levels(
        self,
        lipid_name: str,
        levels: list = None,
        import_rules: dict = default_input_rules,
    ) -> dict:

        if levels is None:
            levels = ["B0"]
        lv_id_dct = {}
        all_lv_id_dct = self.export_all_levels(lipid_name)
        for level in levels:
            if level in supported_levels:
                if level in all_lv_id_dct:
                    lv_id_dct[level] = all_lv_id_dct[level]
                else:
                    raise ValueError(
                        f"Lipid: {lipid_name} cannot be converted into level: {level}. "
                        f"Can be converted into: {all_lv_id_dct}"
                    )
            else:
                raise ValueError(
                    f"Level: {level} not supported. Supported levels: {supported_levels}"
                )

        return lv_id_dct


if __name__ == "__main__":
    t_in_lst = [
        # "GM3(d18:1/18:2(9Z,11Z)(12OH))",
        # "TG P-18:1_18:2(9Z,11Z)(12OH)_18:1(9)(11OH)",
        # "CL(1'-[18:1(9Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",
        # "TG(16:0/18:2/PA)",
        # "PE O-p 32:1",
        # "PE O-a 36:2",
        # "PE O-18:1a/18:1",
        # "PE O-p 36:2",
        # "PE O-18:1p/18:1",
        # "PE-C16:0-C18:2",
        # "HETE",
        # "TG(16:0/18:2/18:2[2xDB,1xOH])",
        # "PC 16:0/18:2[9,12]",
        # "DG 16:2[9,12]_O-18:2_0:0",
        # "14,15-HxB3 (13R)",
        # "C22:5 CE",
        # "15-Keto-PGF2α",
        "PGF2α",
        # "8-iso PGF2a III",
        # "palmitoleic acid",
    ]
    lynx_gen = Encoder(style="COMP_DB")
    for t_in in t_in_lst:
        t1_out = lynx_gen.convert(t_in)
        self.logger.info(f"Input: {t_in} -> Best Output: {t1_out}")
        # t_lv = "B0"
        # t2_out = lynx_gen.export_level(
        #     t_in, level=t_lv, import_rules=default_input_rules
        # )
        # self.logger.info(f"Input: {t_in} -> Output @ Lv {t_lv}: {t2_out}")
        # t_lv_lst = ["B0"]
        # t3_out = lynx_gen.export_levels(
        #     t_in, levels=t_lv_lst, import_rules=default_input_rules
        # )
        # self.logger.info(f"Input: {t_in} -> Output @ Lv {t_lv_lst}: {t3_out}")
        t4_out = lynx_gen.export_all_levels(t_in)
        self.logger.info(f"Input: {t_in} -> Output @ all levels: {t4_out}")

    self.logger.info("fin")
