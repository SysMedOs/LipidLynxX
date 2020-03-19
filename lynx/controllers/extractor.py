# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from typing import Dict, List, Union

from natsort import natsorted
import regex as re

from lynx.models.log import logger
from lynx.models.defaults import (
    rgx_class_dct,
    cv_rgx_dct,
    cv_order_list,
    cv_alias_info,
    default_input_rules,
    default_output_rules,
)
from lynx.controllers.encoder import encode_sub_residues, decode_mod
from lynx.controllers.general_functions import seg_to_str
from lynx.controllers.params_loader import build_input_rules, build_output_rules
from lynx.controllers.parser import rule_parse, parse, parse_mod


class Extractor(object):
    def __init__(self, rules):
        self.rules = rules

    def check_segments(self, lipid_name: str, rule_class: str):
        c = rule_class
        matched_info_dct = {}
        is_this_class = False
        c_search_rgx = self.rules[c].get("SEARCH", None)
        if c_search_rgx.search(lipid_name):
            is_this_class = True
        else:
            if rule_class in ["RESIDUE", "SUM_RESIDUES"]:
                is_this_class = True
            else:
                pass
        c_match_rgx_dct = self.rules[c].get("MATCH", None)
        if is_this_class and isinstance(c_match_rgx_dct, dict):
            for m in c_match_rgx_dct:
                m_pattern = c_match_rgx_dct[m]["MATCH"]
                m_groups = c_match_rgx_dct[m]["GROUPS"]  # type: list
                m_match = m_pattern.match(lipid_name)
                if m_match:
                    matched_dct = {}
                    matched_groups = m_match.capturesdict()
                    matched_info_dct[m] = matched_groups
        return matched_info_dct

    def check_residues(
        self,
        sum_residues: str,
        max_residues: int = 1,
        separator_levels: dict = None,
        separator: str = "-|/",
    ) -> list:

        if separator_levels is None:
            separator_levels = {"B": "", "D": "_", "S": "/"}

        res_lst = re.split(separator, sum_residues)
        res_sep_lst = re.findall(separator, sum_residues)
        if not res_sep_lst:
            res_sep_lst = [""]
        res_sep_levels = []
        for res_sep in res_sep_lst:
            for lv in separator_levels:
                if res_sep == separator_levels[lv]:
                    res_sep_levels.append(lv)

        lv_min = natsorted(res_sep_levels)[0]

        out_res_lst = []
        res_true_lst = []
        if "0:0" in res_lst:
            res_true_lst = [res for res in res_lst if res != "0:0"]

        if len(res_lst) <= max_residues or len(res_true_lst) <= max_residues:
            for res in res_lst:
                # todo: Add rgx here
                matched_info_dct = self.check_segments(res, "RESIDUE")
                for m in matched_info_dct:
                    logger.error(m)
                    matched_dct = matched_info_dct[m]
                    mod_str = matched_dct.get("SUM_MODS", None)
                    # if mod_str:
                    #     mod_info = decode_mod(mod_str)
                    #     logger.info(f"{mod_str}, {mod_info}")
                    logger.error(matched_dct)
                out_res_lst.append(res)

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

        return out_res_lst

    def extract(self, lipid_name: str) -> Dict[str, Union[str, dict]]:

        """
        Main parser to read input abbreviations
        Args:
            lipid_name: input lipid abbreviation to be converted

        Returns:
            parsed_info_dct: parsed information stored as dict

        """

        parsed_info_dct = {}

        for c in self.rules:
            c_lmsd_classes = self.rules[c].get("LMSD_CLASSES", None)
            c_max_res = self.rules[c].get("MAX_RESIDUES", 1)
            res_sep = self.rules[c].get("RESIDUES_SEPARATOR", None)  # type: str
            sep_levels = self.rules[c].get("SEPARATOR_LEVELS", {})  # type: dict
            matched_info_dct = self.check_segments(lipid_name, c)
            for m in matched_info_dct:
                matched_dct = matched_info_dct[m]
                sum_residues_lst = matched_dct.get("SUM_RESIDUES", [])
                if sum_residues_lst:
                    for sum_residues in sum_residues_lst:
                        residues = self.check_residues(
                            sum_residues,
                            max_residues=c_max_res,
                            separator_levels=sep_levels,
                            separator=res_sep,
                        )
                matched_info_dct[m] = {
                    "LMSD_CLASSES": c_lmsd_classes,
                    "SEGMENTS": matched_dct,
                    "RESIDUES_SEPARATOR": res_sep,
                    "SEPARATOR_LEVELS": sep_levels,
                }

                parsed_info_dct[c] = matched_info_dct
            else:
                pass

        if not parsed_info_dct:
            logger.error(f"Failed to decode Lipid: {lipid_name}")

        return parsed_info_dct


if __name__ == "__main__":

    # t_in = "GM3(d18:1/18:2(9Z,12Z))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5Z,8Z,11Z,14Z))"
    t_in = "TG (P-18:1/18:2(9Z,12Z)/5S,15R-DiHETE)"
    extractor = Extractor(rules=default_input_rules)
    t_out = extractor.extract(t_in)

    logger.info(t_out)
    logger.info("FIN")
