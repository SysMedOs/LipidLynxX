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

from typing import Dict, Union, List, Any

from natsort import natsorted
import regex as re

from lynx.controllers.formatter import Formatter
from lynx.models.alias import Alias
from lynx.models.defaults import default_input_rules
from lynx.utils.log import app_logger


class Decoder(object):
    def __init__(self, rules: dict = default_input_rules, logger=app_logger):
        self.rules = rules
        self.formatter = Formatter()
        self.alias = Alias(logger=logger)
        self.logger = logger

    def check_segments(self, lipid_name: str, rule_class: str, rule: str):
        c = rule_class
        matched_info_dct = {}
        is_this_class = False
        c_search_rgx = self.rules[c].get("SEARCH", None)
        if c_search_rgx.search(lipid_name):
            is_this_class = True
        else:
            if rule_class in ["RESIDUE", "SUM_RESIDUES", "RESIDUE_ALIAS", "ALIAS"]:
                is_this_class = True
            else:
                pass
        c_match_rgx_dct = self.rules[c].get("MATCH", None)
        if is_this_class and isinstance(c_match_rgx_dct, dict):
            if rule:
                if rule in c_match_rgx_dct:
                    m_pattern = c_match_rgx_dct[rule]["MATCH"]
                    # m_groups = c_match_rgx_dct[rule]["GROUPS"]  # type: list
                    m_match = m_pattern.match(lipid_name)
                    if m_match:
                        # matched_dct = {}
                        matched_info_dct = m_match.capturesdict()
                else:
                    raise ValueError(f"Can not find rule: {rule} in configuration.")
            else:
                raise ValueError(f"Must provide a {rule} in configuration to search.")

        links = matched_info_dct.get("LINK", None)
        if links:
            if isinstance(links, str):
                if links.lower().startswith(("e", "o", "a")):
                    matched_info_dct["LINK"] = "O-"
                elif links.lower().startswith("p"):
                    matched_info_dct["LINK"] = "P-"
                else:
                    pass
            elif isinstance(links, list):
                std_links = []
                for link in links:
                    if link.lower().startswith(("e", "o", "a")):
                        std_links.append("O-")
                    elif link.lower().startswith("p"):
                        std_links.append("P-")
                matched_info_dct["LINK"] = std_links
            else:
                pass
        else:
            pass
        return matched_info_dct

    def check_alias(self, alias: str, alias_type: str = "RESIDUE") -> str:
        defined_id = ""
        if alias_type.upper().startswith("RESIDUE"):
            lite_alias_info = self.alias.residue_alias
        elif alias_type.upper().startswith("LIPID"):
            lite_alias_info = self.alias.lipid_alias
        else:
            raise ValueError(
                f"Cannot load alias_type {alias_type} from defined_alias.json"
            )

        for alias_rgx in lite_alias_info:
            if re.search(alias_rgx, alias):
                defined_id = lite_alias_info[alias_rgx]
            else:
                pass

        # if not defined_id:
        #     self.logger.debug(
        #         f"Cannot decode alias: {alias} using alias_type: {alias_type}."
        #     )

        return defined_id

    def check_residues(
        self,
        rule: str,
        residues: list,
        sum_residues: str,
        alias=None,
        alias_rule: str = "LipidLynxX.json#LipidLynxX",
        max_residues: int = 1,
        separator_levels: dict = None,
        separator: str = "-|/",
        lmsd_classes: List[str] = None,
    ) -> dict:

        if alias is None:
            alias = []
        if separator_levels is None:
            separator_levels = {"B": "", "D": "_", "S": "/"}

        res_lst = residues
        res_sep_lst = re.findall(separator, sum_residues)
        if not res_sep_lst:
            res_sep_lst = [""]
        for s in res_sep_lst:
            res_lst = [res.strip(s) for res in res_lst]
        res_sep_levels = []
        for res_sep in res_sep_lst:
            for lv in separator_levels:
                if res_sep == separator_levels[lv]:
                    res_sep_levels.append(lv)

        lv_min = natsorted(res_sep_levels)[0]

        out_res_dct = {}
        out_res_lst = []
        res_true_lst = []
        is_fa = False
        for c_lmsd in lmsd_classes:
            if c_lmsd.upper().startswith("FA"):
                is_fa = True
        if not is_fa and "0:0" in res_lst:
            res_true_lst = [res for res in res_lst if res != "0:0"]
        else:
            res_true_lst = res_lst

        if len(res_lst) <= max_residues or len(res_true_lst) <= max_residues:
            for res in res_lst:
                if res in alias:
                    if res != "":
                        def_res = self.check_alias(res, "RESIDUE")
                        if def_res:
                            res = def_res
                    matched_info_dct = self.check_segments(
                        res, "RESIDUE", rule=alias_rule
                    )
                else:
                    matched_info_dct = self.check_segments(res, "RESIDUE", rule=rule)
                # num_o_chk_lst = matched_info_dct.get("NUM_O", [""])
                # if isinstance(num_o_chk_lst, list) and re.match(r'\d?O|O\d?|\d', num_o_chk_lst[0]):
                #     matched_info_dct["MOD_TYPE"] = ["O"] + matched_info_dct.get("MOD_TYPE", [])
                matched_dct = self.formatter.format_residue(matched_info_dct)
                self.logger.debug(f"matched_dct: {matched_dct}")
                out_res_lst.append(res)
                out_res_dct[res] = matched_dct

        if len(out_res_lst) > 1 and lv_min == "B":
            lv_min = "D"

        if lv_min == "D":
            no_res_lst = []
            spb_lst = []
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
            # "RESIDUES_SEPARATOR": res_sep_lst,
            "RESIDUES_SEPARATOR_LEVEL": lv_min,
        }

    def extract_by_class_rule(
        self, lipid_name: str, c: str, lynx_rule_idx: str = "LipidLynxX.json#LipidLynxX"
    ) -> dict:
        c_lmsd_classes = self.rules[c].get("LMSD_CLASSES", None)
        c_max_res = self.rules[c].get("MAX_RESIDUES", 1)
        res_sep = self.rules[c].get("RESIDUES_SEPARATOR", None)  # type: str
        sep_levels = self.rules[c].get("SEPARATOR_LEVELS", {})  # type: dict
        c_rules = self.rules[c].get("MATCH", {})
        matched_info_dct = {}
        for lr in c_rules:
            if re.search(r"Lynx", lr, re.IGNORECASE):
                lynx_rule_idx = lr
        for r in c_rules:
            matched_dct = self.check_segments(lipid_name, c, r)
            sum_residues_lst = matched_dct.get("SUM_RESIDUES", [])
            obs_residues_lst = matched_dct.get("RESIDUE", [])
            alias_lst: list = matched_dct.get("ALIAS", [])
            if sum_residues_lst and len(sum_residues_lst) == 1 and obs_residues_lst:
                residues_dct = self.check_residues(
                    r,
                    obs_residues_lst,
                    sum_residues_lst[0],
                    alias=alias_lst,
                    alias_rule=lynx_rule_idx,
                    max_residues=c_max_res,
                    separator_levels=sep_levels,
                    separator=res_sep,
                    lmsd_classes=c_lmsd_classes,
                )
                # set specific classes into
                for c_lmsd in c_lmsd_classes:
                    if c_lmsd.upper().startswith("FA"):
                        residues_dct["RESIDUES_SEPARATOR_LEVEL"] = "S"
                    elif c_lmsd.upper().startswith("SP"):
                        residues_dct["RESIDUES_SEPARATOR_LEVEL"] = "S"
                        sp_res_abbr_lst = residues_dct["RESIDUES_ORDER"]
                        if len(sp_res_abbr_lst) == 2:
                            if (
                                ";" in sp_res_abbr_lst[1]
                                and ";" not in sp_res_abbr_lst[0]
                            ):
                                residues_dct["RESIDUES_ORDER"] = [
                                    sp_res_abbr_lst[1],
                                    sp_res_abbr_lst[0],
                                ]
                            elif re.match(
                                r"[mdt].*", sp_res_abbr_lst[1]
                            ) and not re.match(r"[mdt].*", sp_res_abbr_lst[0]):
                                residues_dct["RESIDUES_ORDER"] = [
                                    sp_res_abbr_lst[1],
                                    sp_res_abbr_lst[0],
                                ]
                    elif c_lmsd.upper().startswith("ST"):
                        residues_dct["RESIDUES_SEPARATOR_LEVEL"] = "S"
                matched_info_dct[r] = {
                    "LMSD_CLASSES": c_lmsd_classes,
                    "SEGMENTS": matched_dct,
                    "RESIDUES": residues_dct,
                    # "RESIDUES_SEPARATOR": res_sep,
                    # "SEPARATOR_LEVELS": sep_levels,
                }
            elif sum_residues_lst and len(sum_residues_lst) > 1:
                self.logger.error(
                    f"More than two parts of SUM residues matched: {sum_residues_lst}"
                )
            else:
                pass  # nothing found. the rule is not used.
        return matched_info_dct

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
            matched_info_dct = self.extract_by_class_rule(lipid_name, c)
            if matched_info_dct:
                extracted_info_dct[c] = matched_info_dct
            else:
                if lipid_name:
                    def_alias = self.check_alias(lipid_name, "LIPID")
                    if def_alias:
                        self.logger.debug(
                            f"Found Alias: {lipid_name} -> change to {def_alias}"
                        )
                        matched_info_dct = self.extract_by_class_rule(def_alias, c)
                        if matched_info_dct:
                            extracted_info_dct[c] = matched_info_dct

        if not extracted_info_dct:
            self.logger.error(f"Failed to decode Lipid: {lipid_name}")

        return extracted_info_dct


if __name__ == "__main__":

    # LIPID MAPS
    # t_in = "GM3(d18:1/18:2(9Z,12Z))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5Z,8Z,11Z,14Z)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5,8,11,14)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5,8,11,14)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/5S,15R-DiHETE)"

    # MS-DIAL
    # t_in = "TG(16:0/18:2/20:4<OH>)"
    # t_in = "TG(16:0/18:2/HETE)"
    # RefMet
    t_in = "PC(36:2e)"

    extractor = Decoder(rules=default_input_rules)
    t_out = extractor.extract(t_in)

    app_logger.info(t_out)
    app_logger.info("FIN")
