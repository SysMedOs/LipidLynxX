# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List, Union

from natsort import natsorted

from lynx.models.log import logger
from lynx.models.defaults import (
    rgx_class_dct,
    cv_rgx_dct,
    cv_order_list,
    cv_alias_info,
    default_input_rules,
    default_output_rules,
)
from lynx.controllers.encoder import encode_sub_residues
from lynx.controllers.general_functions import seg_to_str
from lynx.controllers.params_loader import build_input_rules, build_output_rules
from lynx.controllers.parser import rule_parse, parse, parse_mod


class Generator(object):
    def __init__(self, export_rules: dict, rule: str):
        self.output_rules = export_rules.get(rule, None)
        self.class_rules = self.output_rules.get("LMSD_CLASSES", {})

    def check_residues(
        self,
        residues: str,
        lmsd_class: str,
        separator_levels: dict = None,
        separator: str = "-|/",
    ) -> str:

        if separator_levels is None:
            separator_levels = {"B": "", "D": "_", "S": "/"}

        res_lst = re.split(separator, residues)
        res_sep_lst = re.findall(separator, residues)
        if not res_sep_lst:
            res_sep_lst = [""]
        res_sep_levels = []
        for res_sep in res_sep_lst:
            for lv in separator_levels:
                if res_sep == separator_levels[lv]:
                    res_sep_levels.append(lv)

        lv_min = natsorted(res_sep_levels)[0]

        lmsd_rules = self.output_rules.get("LMSD_CLASSES", {}).get(lmsd_class, {})
        c_max_res_count = lmsd_rules.get("MAX_RESIDUES", None)
        sep_rules = self.output_rules.get("SEPARATORS")
        out_sep = sep_rules.get("SEPARATOR_LEVELS", {}).get(lv_min, "")

        out_res_lst = []
        res_true_lst = []
        if "0:0" in res_sep_lst:
            res_true_lst = [res for res in res_lst if res != "0:0"]
        if len(res_lst) <= c_max_res_count or len(res_true_lst) <= c_max_res_count:
            for res in res_lst:
                res_info = parse(res)
                res_lynx_code = encode_sub_residues(res)
                logger.info(res_lynx_code)
                logger.info(res_info)
                out_res_lst.append(res_lynx_code)

        if lv_min == "D":
            o_lst = []
            p_lst = []
            r_lst = []
            for r in out_res_lst:
                if r.startswith("O-"):
                    o_lst.append(r)
                elif r.startswith("P-"):
                    p_lst.append(r)
                else:
                    r_lst.append(r)

            out_res_lst = natsorted(o_lst) + natsorted(p_lst) + natsorted(r_lst)

        return out_sep.join(out_res_lst)

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
                    out_seg_lst.append(patterns_dct[s_rgx])
            out_seg_lst = list(filter(None, list(set(out_seg_lst))))
        if out_seg_lst:
            out_seg_str = out_seg_lst[0]
        else:
            out_seg_str = ""
        return out_seg_str

    def check_segments(self, parsed_info: dict, input_rule: str):
        segments_dct = {}
        lmsd_classes = parsed_info.get("LMSD_CLASSES", None)
        segments = parsed_info["SEGMENTS"]
        res_sep_str = parsed_info["RESIDUES_SEPARATOR"]
        res_sep_levels = parsed_info["SEPARATOR_LEVELS"]
        c_str = segments.get("CLASS", "")
        for c in lmsd_classes:
            if c in self.class_rules:
                c_dct = segments_dct.get(c, {})
                lc_class_str = self.check_rest(segments.get("CLASS", ""), "CLASS", c)
                if lc_class_str:
                    out_sum_res_str = self.check_residues(
                        segments.get("SUM_RESIDUES", ""),
                        lmsd_class=c,
                        separator_levels=res_sep_levels,
                        separator=res_sep_str,
                    )
                    c_dct[input_rule] = {
                        "LMSD_CLASS": c,
                        "SEGMENTS": {
                            "PREFIX": self.check_rest(
                                segments.get("PREFIX", ""), "PREFIX", c
                            ),
                            "CLASS": lc_class_str,
                            "SUFFIX": self.check_rest(
                                segments.get("SUFFIX", ""), "SUFFIX", c
                            ),
                            "BRACKET_LEFT": segments.get("BRACKET_LEFT", "("),
                            "SUM_RESIDUES": out_sum_res_str,
                            "BRACKET_RIGHT": segments.get("BRACKET_RIGHT", ")"),
                        },
                    }
                    segments_dct[c] = c_dct
                    logger.info(
                        f'Assign lipid class {lc_class_str} to {segments.get("CLASS", "")} from class {c}'
                    )

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
        if len(candidates_lst) > 1:
            code_len = 1
            for tmp_code in candidates_lst:
                tmp_len = len(tmp_code)
                if max(code_len, tmp_len) == tmp_len:
                    best_abbr = tmp_code
        elif len(candidates_lst) == 1:
            best_abbr = candidates_lst[0]
        else:
            logger.warning("Failed to generate abbreviation for this lipid...")

        return best_abbr

    def export(self, lipid_name: str, import_rules: dict = default_input_rules):

        parsed_info = rule_parse(lipid_name, rules=import_rules)
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
    import os

    from lynx.controllers.general_functions import load_folder
    from lynx.models.rules import InputRules

    check_in_rules_folder = r"../configurations/rules/input"
    file_path_lst = load_folder(check_in_rules_folder, file_type=".json")
    logger.debug(f"Fund JSON config files: \n {file_path_lst}")

    in_rule_tests = {}
    for f in file_path_lst:
        tmp_rules = InputRules(f)
        idx_lst = [os.path.basename(f)] + tmp_rules.source
        idx = "#".join(idx_lst)
        rules_info = tmp_rules.raw_rules.get("LIPID_CLASSES", {})
        tmp_tests = {}
        for c in rules_info:
            in_rule_tests[f"{c}@{idx}"] = rules_info[c].get("EXAMPLES", [])

    lynx_gen = Generator(export_rules=default_output_rules, rule="LipidLynxX@20200214")
    test_sum = {}
    for i in in_rule_tests:
        i_lst = in_rule_tests[i]
        i_test_dct = {}
        for e in i_lst:
            o = lynx_gen.export(e, import_rules=default_input_rules)
            i_test_dct[e] = o
        test_sum[i] = i_test_dct
    for k in test_sum:
        logger.info(f"Results of {k}:")
        for r in test_sum[k]:
            r_out = test_sum[k].get(r, "")
            if len(r_out) >= len(r):
                logger.info(f"Input: {r} -> Output: {r_out}")
            elif r_out.endswith(")") or r_out.endswith(">"):
                logger.info(f"Input: {r} -> Output: {r_out}")
            else:
                logger.warning(f"{k} Input: {r} -> Output: {r_out}")

    logger.info("fin")
