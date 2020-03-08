# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Dict, List, Union

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

    def check_lmsd_class(self, parsed_info: dict, input_rule: str):
        lmsd_info = {}
        lmsd_classes = parsed_info.get("LMSD_CLASSES", None)
        segments = parsed_info["SEGMENTS"]
        res_sep_str = parsed_info["RESIDUES_SEPARATOR"]
        res_sep_levels = parsed_info["SEPARATOR_LEVELS"]
        c_str = segments.get("CLASS", "")
        for c in lmsd_classes:
            if c in self.class_rules:
                lmsd_patterns_lst = self.class_rules[c].get("INPUT_PATTERNS")
                for c_rgx in lmsd_patterns_lst:
                    logger.debug(f"Test {c_str} on LMSD: {c} using {c_rgx}")
                    c_searched = c_rgx.search(c_str)
                    if c_searched:
                        c_dct = lmsd_info.get(c, {})
                        c_dct[input_rule] = {
                            "SEGMENTS": segments,
                            "RESIDUES_SEPARATOR": res_sep_str,
                        }
                        temp_out_str = self.check_residues(
                            segments,
                            lmsd_class=c,
                            separator=res_sep_str,
                        )
        return lmsd_info

    def check_residues(
            self, segments: dict, lmsd_class: str, separator: str = "-|/"
    ) -> str:

        output_str = ""
        sum_res_str = segments.get("SUM_RESIDUES", "")
        lc_class_str = segments.get("CLASS", "")
        res_lst = re.split(separator, sum_res_str)
        res_sep_lst = re.findall(separator, sum_res_str)
        logger.info(res_sep_lst)
        lmsd_rules = self.output_rules.get("LMSD_CLASSES", {}).get(lmsd_class, {})
        sep_rules = self.output_rules.get("SEPARATORS")
        c_max_res_count = lmsd_rules.get("MAX_RESIDUES", None)
        if len(res_sep_lst) <= c_max_res_count:
            for res in res_lst:
                res_info = parse(res)
                res_lynx_code = encode_sub_residues(res)
                logger.info(res_lynx_code)
                logger.info(res_info)
        return output_str

    def export(self, lipid_name: str, import_rules: dict = default_input_rules):

        parsed_info = rule_parse(lipid_name, rules=import_rules)
        for p in parsed_info:
            p_info = parsed_info[p]
            logger.info(p_info)
            for in_r in p_info:
                r_info = p_info[in_r]  # type: dict
                lmsd_info = self.check_lmsd_class(r_info, in_r)





if __name__ == "__main__":

    # examples = [
    #     "LPE 16:0",
    #     "LPE 16:0/0:0",
    #     "LPE O-18:1",
    #     "LPE P-16:0",
    #     "PE 34:2",
    #     "PE 16:0_18:2",
    #     "PE O-34:2",
    #     "PE O-16:0_18:2",
    #     "PE P-34:2",
    #     "PE P-16:0_18:2",
    #     "PIP 34:1",
    #     "PIP 16:0_18:1",
    #     "LPIP 16:0",
    #     "LPIP 16:0/0:0",
    #     "PIP2 34:2",
    #     "PIP2 16:0_18:2",
    #     "LPIP2 16:0",
    #     "LPIP2 16:0/0:0",
    #     "PIP3 34:2",
    #     "PIP3 16:0_18:2",
    #     "LPIP3 16:0",
    #     "LPIP3 16:0/0:0",
    #     "PEtOH 34:2",
    #     "PEtOH 16:0_18:2",
    #     "BMP 34:1",
    #     "BMP 16:0_18:1",
    #     "MG(16:0)",
    #     "MG(0:0/0:0/16:0)",
    #     "MG(0:0/16:0/0:0)",
    #     "MG(16:0/0:0/0:0)",
    #     "DG(34:2)",
    #     "DG(16:0_18:2)",
    #     "DG(O-34:2)",
    #     "DG(O-16:0_18:2)",
    #     "DG(P-34:2)",
    #     "DG(P-16:0_18:2)",
    #     "DG(P-16:0/0:0/18:2)",
    #     "TG(52:2)",
    #     "TG(16:0_18:0_18:2)",
    #     "TG(16:0/18:2/18:0)",
    #     "TG(O-52:2)",
    #     "TG(O-16:0_18:0_18:2)",
    #     "TG(O-16:0/18:2/18:0)",
    #     "TG(P-52:2)",
    #     "TG(P-16:0_18:0_18:2)",
    #     "TG(P-16:0/18:2/18:0)",
    # ]

    examples = ["TG(P-16:0_18:2(9Z,11Z)_18:1(9Z))"]

    for e in examples:
        o = generate(
            e,
            import_rules=default_input_rules,
            export_rules=default_output_rules,
            rule="LipidLynxX@20200214",
        )
        logger.debug(e)
        logger.debug(o)

    logger.info("fin")
