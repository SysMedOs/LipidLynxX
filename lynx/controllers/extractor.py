# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from typing import Dict, Union

from natsort import natsorted
import regex as re

from lynx.controllers.formatter import Formatter
from lynx.models.defaults import default_input_rules
from lynx.utils.log import logger


class Extractor(object):
    def __init__(self, rules: dict = default_input_rules):
        self.rules = rules
        self.formatter = Formatter()

    def check_segments(self, lipid_name: str, rule_class: str, rule: str):
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
            if rule:
                if rule in c_match_rgx_dct:
                    m_pattern = c_match_rgx_dct[rule]["MATCH"]
                    m_groups = c_match_rgx_dct[rule]["GROUPS"]  # type: list
                    m_match = m_pattern.match(lipid_name)
                    if m_match:
                        matched_dct = {}
                        matched_groups = m_match.capturesdict()
                        matched_info_dct = matched_groups
                else:
                    raise ValueError(f"Can not find rule: {rule} in configuration.")
            else:
                raise ValueError(f"Must provide a {rule} in configuration to search.")
                # for m in c_match_rgx_dct:
                #     m_pattern = c_match_rgx_dct[m]["MATCH"]
                #     m_groups = c_match_rgx_dct[m]["GROUPS"]  # type: list
                #     m_match = m_pattern.match(lipid_name)
                #     if m_match:
                #         matched_dct = {}
                #         matched_groups = m_match.capturesdict()
                #         matched_info_dct[m] = matched_groups
        return matched_info_dct

    def check_residues(
        self,
        rule: str,
        sum_residues: str,
        max_residues: int = 1,
        separator_levels: dict = None,
        separator: str = "-|/",
    ) -> dict:

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

        out_res_dct = {}
        out_res_lst = []
        res_true_lst = []
        if "0:0" in res_lst:
            res_true_lst = [res for res in res_lst if res != "0:0"]

        if len(res_lst) <= max_residues or len(res_true_lst) <= max_residues:
            for res in res_lst:
                # todo: Add rgx here
                matched_info_dct = self.check_segments(res, "RESIDUE", rule=rule)
                matched_dct = self.formatter.format_residue(matched_info_dct)
                logger.error(matched_dct)
                out_res_lst.append(res)
                out_res_dct[res] = matched_dct

        if lv_min == "D":
            no_res_lst = []
            o_lst = []
            p_lst = []
            r_lst = []
            for r in out_res_lst:
                if r == "0:0":
                    no_res_lst.append(r)
                elif r.startswith("O-"):
                    o_lst.append(r)
                elif r.startswith("P-"):
                    p_lst.append(r)
                else:
                    r_lst.append(r)

            out_res_lst = (
                natsorted(no_res_lst)
                + natsorted(o_lst)
                + natsorted(p_lst)
                + natsorted(r_lst)
            )

        return {
            "RESIDUES_ORDER": out_res_lst,
            "RESIDUES_INFO": out_res_dct,
            "RESIDUES_SEPARATOR": res_sep_lst,
            "RESIDUES_SEPARATOR_LEVEL": lv_min,
        }

    def extract(self, lipid_name: str) -> Dict[str, Union[str, dict]]:

        """
        Main parser to read input abbreviations
        Args:
            lipid_name: input lipid abbreviation to be converted

        Returns:
            extracted_info_dct: parsed information stored as dict

        """

        extracted_info_dct = {}

        for c in self.rules:
            c_lmsd_classes = self.rules[c].get("LMSD_CLASSES", None)
            c_max_res = self.rules[c].get("MAX_RESIDUES", 1)
            res_sep = self.rules[c].get("RESIDUES_SEPARATOR", None)  # type: str
            sep_levels = self.rules[c].get("SEPARATOR_LEVELS", {})  # type: dict

            c_rules = self.rules[c].get("MATCH", {})
            matched_info_dct = {}
            for r in c_rules:
                matched_dct = self.check_segments(lipid_name, c, r)
                sum_residues_lst = matched_dct.get("SUM_RESIDUES", [])
                if sum_residues_lst and len(sum_residues_lst) == 1:
                    residues_dct = self.check_residues(
                        r,
                        sum_residues_lst[0],
                        max_residues=c_max_res,
                        separator_levels=sep_levels,
                        separator=res_sep,
                    )
                    matched_info_dct[r] = {
                        "LMSD_CLASSES": c_lmsd_classes,
                        "SEGMENTS": matched_dct,
                        "Residues": residues_dct,
                        # "RESIDUES_SEPARATOR": res_sep,
                        # "SEPARATOR_LEVELS": sep_levels,
                    }
                elif sum_residues_lst and len(sum_residues_lst) > 1:
                    raise ValueError(
                        f"More than two parts of SUM residues matched: {sum_residues_lst}"
                    )
                else:
                    pass  # nothing found. the rule is not used.

            if matched_info_dct:
                extracted_info_dct[c] = matched_info_dct

        if not extracted_info_dct:
            logger.error(f"Failed to decode Lipid: {lipid_name}")

        return extracted_info_dct


if __name__ == "__main__":

    # t_in = "GM3(d18:1/18:2(9Z,12Z))"
    t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5Z,8Z,11Z,14Z)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5,8,11,14)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5,8,11,14)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/5S,15R-DiHETE)"
    extractor = Extractor(rules=default_input_rules)
    t_out = extractor.extract(t_in)

    logger.info(t_out)
    logger.info("FIN")
