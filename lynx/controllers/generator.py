# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import List

from natsort import natsorted

from lynx.controllers.encoder import encode_sub_residues
from lynx.controllers.parser import rule_parse, parse
from lynx.controllers.params_loader import load_output_rule
from lynx.controllers.extractor import Extractor
from lynx.models.residue import Residue
from lynx.models.defaults import default_output_rules, default_input_rules
from lynx.utils.log import logger


class Generator(object):
    def __init__(self, output_rules: dict = default_output_rules, rule: str = "LipidLynxX",
                 input_rules: dict = default_input_rules):
        self.output_rules = load_output_rule(output_rules, rule)
        self.class_rules = self.output_rules.get("LMSD_CLASSES", {})
        self.mod_rules = self.output_rules.get("MODS", {}).get("MOD", {})
        self.sum_mod_rules = self.output_rules.get("MODS", {}).get("SUM_MODS", {})
        self.residue_rules = self.output_rules.get("RESIDUES", {}).get("RESIDUE", {})
        self.separator_levels = self.output_rules.get("SEPARATORS", {}).get("SEPARATOR_LEVELS", {})
        self.extractor = Extractor(rules=input_rules)

    @staticmethod
    def get_best_candidate(candidates: List[str]) -> str:
        best_candidate = ""
        if len(candidates) > 1:
            code_len = 1
            for tmp_code in candidates:
                tmp_len = len(tmp_code)
                if max(code_len, tmp_len) == tmp_len:
                    best_candidate = tmp_code
                    code_len = tmp_len
        elif len(candidates) == 1:
            best_candidate = candidates[0]
        else:
            logger.warning("Failed to generate abbreviation for this lipid...")
        return best_candidate

    def check_rest(self, segment_text: str, segment_name: str, lmsd_class: str):
        patterns_dct = self.class_rules[lmsd_class].get(segment_name)
        out_seg_lst = []
        if segment_text and patterns_dct:
            for s_rgx in patterns_dct:
                logger.debug(
                    f"Test {segment_text} on {segment_name} of {lmsd_class} using {s_rgx}"
                )
                s_matched = s_rgx.match(segment_text)
                if s_matched:
                    defined_seg = patterns_dct[s_rgx]
                    if defined_seg == "EXCEPTIONS":
                        if lmsd_class in ["GP12"]:
                            out_seg_lst.append(segment_text)
                        if lmsd_class in ["SP05", "SP06"]:
                            out_seg_lst.append(segment_text)
                    else:
                        out_seg_lst.append(defined_seg)

            out_seg_lst = list(filter(None, list(set(out_seg_lst))))

        return self.get_best_candidate(out_seg_lst)

    def get_residues(self, residues: dict):
        residues_order = residues.get("RESIDUES_ORDER", [])
        residues_sep_level = residues.get("RESIDUES_SEPARATOR_LEVEL", "")
        residues_info = residues.get("RESIDUES_INFO", [])
        res_count = len(list(residues_info.keys()))
        sum_residues_str = ''
        res_lv_id_dct = {}
        res_lv_dct = {}
        sum_lv_lst = []
        for res_abbr in residues_info:
            res_obj = Residue(residues_info[res_abbr])
            res_lv_id_dct[res_abbr] = res_obj.linked_ids
            res_lv_dct[res_abbr] = list(res_obj.linked_ids.keys())
            sum_lv_lst.extend(res_lv_dct[res_abbr])

        sum_lv_lst = natsorted(set(sum_lv_lst))
        for res in res_lv_id_dct:
            r_lv_dct = res_lv_id_dct[res]
            if list(r_lv_dct.keys()) == ['0']:
                for lv in sum_lv_lst:
                    r_lv_dct[lv] = r_lv_dct["0"]
        sum_lv_id_dct = {}
        for sum_lv in sum_lv_lst:
            lv_id_lst = []
            for r in residues_order:
                r_lv_id = res_lv_id_dct.get(r, {}).get(sum_lv, None)
                if r_lv_id:
                    lv_id_lst.append(r_lv_id)
            if len(lv_id_lst) == res_count:
                sum_lv_id_dct[sum_lv] = lv_id_lst

        sum_res_sep_lv_lst = []
        if residues_sep_level == 'S':
            sum_res_sep_lv_lst = ['D', 'S']
        elif residues_sep_level == 'D':
            sum_res_sep_lv_lst = ['D']
        elif residues_sep_level == 'B':
            sum_res_sep_lv_lst = ['B']

        sum_res_id_lv_dct = {}

        for sep_lv in sum_res_sep_lv_lst:
            for res_lv in sum_lv_id_dct:
                sum_res_id_lv_dct[f'{sep_lv}{res_lv}'] = f'{self.separator_levels[sep_lv]}'.join(sum_lv_id_dct[res_lv])

        return sum_res_id_lv_dct

    def check_segments(self, parsed_info: dict, input_rule: str):
        segments_dct = {}
        lmsd_classes = parsed_info.get("LMSD_CLASSES", None)
        segments = parsed_info["SEGMENTS"]
        residues = parsed_info.get("RESIDUES", {})
        sum_residues_str = self.get_residues(residues)
        for c in lmsd_classes:
            if c in self.class_rules:
                pass

        return segments_dct

    def compile_segments(self, segments: dict):
        comp_seg_dct = {}
        for c in segments:
            c_seg_dct = segments[c]
            c_comp_seg_dct = comp_seg_dct.get(c, {})
            for input_rule in c_seg_dct:
                c_seg = c_seg_dct[input_rule].get("SEGMENTS", {})
                c_out_orders = self.output_rules["LMSD_CLASSES"][c].get("ORDER", [])
                ordered_seg_lst = []
                for order in c_out_orders:
                    if order in c_seg:
                        if c_seg[order]:
                            ordered_seg_lst.append(c_seg[order])
                        else:
                            pass
                    else:
                        logger.debug(f"Cannot find segment: {order} in {c_seg}")
                c_comp_seg_dct[input_rule] = "".join(ordered_seg_lst)
            comp_seg_dct[c] = c_comp_seg_dct

        return comp_seg_dct

    def get_best_abbreviation(self, candidates_lst: List[dict]) -> str:
        best_abbr = ""
        chk_abbr_dct = {}
        chk_abbr_lst = []
        for c_info in candidates_lst:
            for c in c_info:
                c_abbr_lst = chk_abbr_dct.get(c, [])
                for in_rule in c_info[c]:
                    chk_abbr_lst.append(c_info[c][in_rule])
                    c_abbr_lst.append(c_info[c][in_rule])
                    chk_abbr_dct[c] = c_abbr_lst
        candidates_lst = list(filter(None, list(set(chk_abbr_lst))))
        for abbr in candidates_lst:
            if abbr.startswith("FA"):
                if re.search(r":\d[oep]", abbr) or re.search(
                        r"[OP]-\d{1,2}:\da?", abbr
                ):
                    candidates_lst.remove(abbr)
                    candidates_lst.append(abbr[2:])
        candidates_lst = list(filter(None, list(set(chk_abbr_lst))))

        return self.get_best_candidate(candidates_lst)

    def export(self, lipid_name: str, import_rules: dict = default_input_rules):

        parsed_info = self.extractor.extract(lipid_name)
        export_info = []
        for p in parsed_info:
            p_info = parsed_info[p]
            logger.info(p_info)
            for in_r in p_info:
                r_info = p_info[in_r]  # type: dict
                lmsd_info = self.check_segments(r_info, in_r)
                comp_dct = self.compile_segments(lmsd_info)
                export_info.append(comp_dct)
        best_export_str = self.get_best_abbreviation(export_info)
        logger.debug(export_info)
        logger.debug(best_export_str)

        return best_export_str


if __name__ == "__main__":
    t_in = "GM3(d18:1/18:2(9Z,11Z)(12OH))"
    lynx_gen = Generator()
    t_out = lynx_gen.export(t_in, import_rules=default_input_rules)
    logger.warning(f"Input: {t_in} -> Output: {t_out}")

    logger.info("fin")
