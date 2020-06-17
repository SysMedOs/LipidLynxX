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

from typing import Union

import regex as re

from lynx.utils.log import app_logger
from lynx.utils.file_handler import get_json


class InputRules(object):
    """
    Read input rules from json file and build corresponding regular expressions
    """

    def __init__(self, rules: Union[str, dict], logger=app_logger):
        if isinstance(rules, dict):
            pass
        elif isinstance(rules, str):
            rules = get_json(rules)
        else:
            raise TypeError
        self.raw_rules = rules
        self.sources = self.raw_rules["SOURCES"]
        self.date = self.raw_rules.get("_DATE", 20200214)
        self.authors = self.raw_rules.get("_AUTHORS", ["example@uni-example.de"])
        mods = self.raw_rules.get("MODS", {})
        if mods:
            self.supported_mods = list(mods.keys())
        else:
            self.supported_mods = []
        self.supported_residues = list(self.raw_rules["RESIDUES"].keys())
        self.supported_classes = list(self.raw_rules["LIPID_CLASSES"].keys())
        self.supported_keys = (
            self.supported_mods + self.supported_residues + self.supported_classes
        )
        self.separators = self.raw_rules["SEPARATORS"]
        self.rules = {}
        self.is_structure_valid = self.__check__()

        if self.is_structure_valid:
            pass
        else:
            raise ValueError(
                "Cannot find valid sections for key [SOURCE, RESIDUES, LIPID_CLASSES]"
            )
        self.rules = self.build()
        self.is_validated = self.validate()
        self.logger = logger
        # self.logger.info(
        #     f"Load input rule: {self.sources}\n"
        #     f"Last modified: {self.date}\n"
        #     f"Authors: {self.authors}"
        # )

    def __replace_refs__(self, rules: dict):
        ref_rgx = re.compile(r"(\$)((\.[A-Z_]{1,99})+)(\.0)")
        for rule in rules:
            rule_dct = rules[rule]
            pattern = rule_dct["PATTERN"]
            # self.logger.debug(f"Pattern: {pattern}")
            ref_lst = re.findall(ref_rgx, pattern)
            ref_replace_dct = {}
            if ref_lst:
                # self.logger.info(f"Found Refs: {ref_lst}")
                for ref_tpl in ref_lst:
                    ref_seg_dct = {}
                    if len(ref_tpl) >= 2:
                        ref_seg_lst = ref_tpl[1].split(".")
                        ref_seg_lst = [s for s in ref_seg_lst if s != "" and s != " "]
                        if ref_seg_lst:
                            ref_pattern_lst = []
                            ref_info = ""
                            for ref in ref_seg_lst:
                                if not ref_seg_dct:
                                    if ref in rules:
                                        ref_info = rules[ref]["PATTERN"]
                                        if ref_info and isinstance(ref_info, str):
                                            # self.logger.info(f"Found Ref Pattern: {ref_lst}")
                                            ref_patt = (
                                                r"\$\."
                                                + "\\.".join(ref_seg_lst)
                                                + r"\.0"
                                            )
                                            ref_replace_dct[ref_patt] = ref_info
                                            break
                                    else:
                                        if ref in self.raw_rules:
                                            ref_seg_dct = self.raw_rules[ref]
                                            ref_pattern_lst.append(r"\." + ref)
                                else:
                                    if ref in rules:
                                        ref_info = rules[ref]["PATTERN"]
                                        if ref_info and isinstance(ref_info, str):
                                            # self.logger.info(f"Found Ref Pattern: {ref_lst}")
                                            ref_patt = (
                                                r"\$\." + "\\.".join(ref_seg_lst) + ".0"
                                            )
                                            ref_replace_dct[ref_patt] = ref_info
                                            break
                                    elif ref in ref_seg_dct:
                                        ref_info = ref_seg_dct[ref]
                                        if isinstance(ref_info, str):
                                            ref_pattern_lst.append(r"\." + ref)
                                            ref_pattern = (
                                                r"\$" + "".join(ref_pattern_lst) + ".0"
                                            )
                                            ref_replace_dct[ref_pattern] = ref_info
                                            # self.logger.info(f"Found Ref Pattern: {ref_lst}")
                                            break
                                        elif isinstance(ref_info, dict):
                                            ref_seg_dct = ref_info
                                        else:
                                            raise KeyError
                                    else:
                                        raise KeyError
                        else:
                            raise ValueError
                    else:
                        raise ValueError
                # self.logger.info(f"Replace Refs: {ref_replace_dct}")
                replace_match = False
                for ref_replace in ref_replace_dct:
                    replaced_pattern = re.sub(
                        ref_replace,
                        str(ref_replace_dct[ref_replace]),
                        rule_dct["PATTERN"],
                    )
                    rule_dct["PATTERN"] = replaced_pattern
                    # self.logger.info(
                    #     f"Replaced pattern to {replaced_pattern} by {ref_replace}"
                    # )
                    replace_match = True
                if replace_match:
                    rule_dct["MATCH"] = re.compile(rule_dct["PATTERN"])
        return rules

    def __check__(self):
        is_structure_valid = False
        if (
            isinstance(self.sources, list)
            and isinstance(self.supported_residues, list)
            and isinstance(self.supported_classes, list)
        ):
            is_structure_valid = True
        else:
            pass

        return is_structure_valid

    def __build__(self, lipid_classes: list, group: str) -> dict:
        """
        Build regular expressions from segments of rules

        Args:
            lipid_classes: List of supported lipid classes or residues
            group: define the rules is for "Lipid_CLASSES" or "RESIDUES"

        Returns:
            a dict of compiled regular expression patterns

        """

        rules = {}
        for c in lipid_classes:
            temp_c_dct = self.raw_rules.get(group, {}).get(c, {})
            if temp_c_dct:
                order_lst = temp_c_dct.get("ORDER", [])
                optional_lst = temp_c_dct.get("OPTIONAL", [])
                repeat_dct = temp_c_dct.get("REPEAT", {})
                pattern_str = ""
                seg_lst = []
                for seg in order_lst:
                    is_optional = False
                    if seg in optional_lst:
                        is_optional = True
                    seg_str = temp_c_dct.get(seg, "")
                    if c in self.supported_keys and seg_str in self.supported_residues:
                        seg_str = self.rules[seg_str].get("PATTERN", "")
                        if seg_str:
                            if seg in repeat_dct:
                                is_optional = False
                                repeat_info = repeat_dct[seg]
                                if (
                                    isinstance(repeat_info, list)
                                    and len(repeat_info) == 2
                                ):
                                    pattern_str += (
                                        f"({seg_str})"
                                        + "{"
                                        + f"{repeat_info[0]},{repeat_info[1]}"
                                        + "}"
                                    )
                                else:
                                    pattern_str += f"({seg_str})*"
                            else:
                                pattern_str += seg_str
                            seg_lst.append(seg)
                        else:
                            raise KeyError(f"Rule pattern {seg_str} not defined!")
                    else:
                        if seg in repeat_dct:
                            is_optional = False
                            repeat_info = repeat_dct[seg]
                            if isinstance(repeat_info, list) and len(repeat_info) == 2:
                                pattern_str += (
                                    f"(?P<{seg}>{seg_str})"
                                    + "{"
                                    + f"{repeat_info[0]},{repeat_info[1]}"
                                    + "}"
                                )
                            else:
                                if seg in optional_lst:
                                    pattern_str += f"(?P<{seg}>{seg_str})*"
                                else:
                                    pattern_str += f"(?P<{seg}>{seg_str})+"
                        else:
                            pattern_str += f"(?P<{seg}>{seg_str})"
                        seg_lst.append(seg)

                    if is_optional:
                        pattern_str += "?"

                rules[c] = {
                    "GROUPS": seg_lst,
                    "PATTERN": pattern_str,
                    "MATCH": re.compile(pattern_str),
                    "CLASS": temp_c_dct.get("CLASS", c),
                    "LMSD_CLASSES": temp_c_dct.get("LMSD_CLASSES", [c]),
                    "MAX_RESIDUES": temp_c_dct.get("MAX_RESIDUES", 1),
                    "RESIDUES_SEPARATOR": self.separators.get(
                        "RESIDUES_SEPARATOR", "_|/"
                    ),
                }
        return rules

    def build(self) -> dict:
        sum_rules = self.__build__(self.supported_mods, "MODS")
        sum_rules.update(self.__build__(self.supported_residues, "RESIDUES"))
        # sum_rules = self.__replace_refs__(sum_rules)
        sum_rules.update(self.__build__(self.supported_classes, "LIPID_CLASSES"))
        sum_rules = self.__replace_refs__(sum_rules)
        self.rules = sum_rules
        return sum_rules

    def validate(self) -> bool:
        """
        Validate if the regular expression fits to all examples in the JSON config file.

        Returns:
            True if all test passed, otherwise return False.

        """
        if not self.rules:
            self.rules = self.build()
        else:
            pass
        valid_lst = []
        valid_dct = {}
        for c in self.supported_keys:

            if c in self.supported_mods:
                temp_c_dct = self.raw_rules.get("MODS", {}).get(c, {})
                seg_typ = "MODS"
            elif c in self.supported_residues:
                temp_c_dct = self.raw_rules.get("RESIDUES", {}).get(c, {})
                seg_typ = "RESIDUES"
            elif c in self.supported_classes:
                temp_c_dct = self.raw_rules.get("LIPID_CLASSES", {}).get(c, {})
                seg_typ = "LIPID_CLASSES"
            else:
                temp_c_dct = self.raw_rules.get("LIPID_CLASSES", {}).get(c, {})
                seg_typ = "LIPID_CLASSES"
            pattern_dct = self.rules.get(c, "")
            if isinstance(pattern_dct, dict):
                pattern_str = pattern_dct.get("PATTERN", None)
            else:
                pattern_str = None
            if pattern_str and temp_c_dct:

                max_res_count = self.raw_rules[seg_typ][c].get("MAX_RESIDUES", 1)
                # self.logger.debug(f"Validating {c} pattern: {pattern_str}")
                c_pattern = re.compile(pattern_str)
                test_lst = temp_c_dct.get("EXAMPLE", [])
                test_dct = {}
                for test_str in test_lst:
                    fit_class = False
                    class_rgx = re.compile(pattern_dct["CLASS"])
                    class_search = class_rgx.search(test_str)
                    # if class_search:
                    #     self.logger.debug(f"{test_str} fit to class {c}")
                    # else:
                    #     self.logger.error(f"{test_str} NOT fit to class {c}")
                    c_match = c_pattern.match(test_str)
                    if c_match:
                        c_matched_res_dct = c_match.groupdict()
                        # self.logger.debug(c_matched_res_dct)
                        if "SUM_RESIDUES" in c_matched_res_dct:
                            c_sum_res = c_matched_res_dct["SUM_RESIDUES"]
                            c_sum_res_lst = re.split(
                                self.separators.get("RESIDUES_SEPARATOR", "_|/"),
                                c_sum_res,
                            )
                            res_pat_str = self.rules.get("RESIDUE", {}).get(
                                "PATTERN", None
                            )
                            res_pat = re.compile(res_pat_str)
                            if 0 < len(c_sum_res_lst) <= max_res_count and res_pat_str:
                                for res in c_sum_res_lst:
                                    if res_pat.match(res):
                                        valid_lst.append(True)
                                        test_dct[f"{test_str} # {res}"] = "passed"
                                        self.logger.debug(f"Check: {test_str} # {res}")
                                        self.logger.debug(res_pat.match(res).groupdict())
                                    else:
                                        valid_lst.append(False)
                                        test_dct[f"{test_str} # {res}"] = "failed"
                                        self.logger.error(f"{test_str} # {res} failed")
                            else:
                                raise ValueError(
                                    f"Number of residues exceed the defined limit: "
                                    f"{max_res_count} for Lipid_Class: {c} from source: {self.sources}"
                                )
                        else:
                            valid_lst.append(True)
                            test_dct[test_str] = "passed"
                    else:
                        valid_lst.append(False)
                        test_dct[test_str] = "failed"
                        self.logger.error(f"{test_str} failed")

                valid_dct[c] = test_dct

        if False in valid_lst:
            is_valid = False
        else:
            is_valid = True

        return is_valid


class OutputRules(object):
    """
    Read output rules from json file and build corresponding regular expressions
    """

    def __init__(self, rules: Union[str, dict]):
        if isinstance(rules, dict):
            pass
        elif isinstance(rules, str):
            rules = get_json(rules)
        else:
            raise TypeError
        self.raw_rules = rules
        self.date = self.raw_rules.get("_DATE", 20200214)
        self.authors = self.raw_rules.get("_AUTHORS", ["example@uni-example.de"])
        self.nomenclature = self.raw_rules.get("NOMENCLATURE", "LipidLynxX")
        self.supported_lmsd_classes = list(self.raw_rules["LMSD_CLASSES"].keys())
        self.separators = self.raw_rules["SEPARATORS"]
        self.mods = self.raw_rules.get("MODS", {})
        self.residues = self.raw_rules.get("RESIDUES", {})
        self.rules = self.build()
        self.is_structure_valid = self.__check__()
        # self.logger.info(
        #     f"Load rules for nomenclature {self.nomenclature}\n"
        #     f"Last modified: {self.date}\n"
        #     f"Authors: {self.authors}"
        # )

    def __check__(self):
        is_structure_valid = False
        if isinstance(self.nomenclature, str) and isinstance(
            self.supported_lmsd_classes, list
        ):
            is_structure_valid = True
        else:
            pass

        return is_structure_valid

    def build(self):
        rules = {
            "LMSD_CLASSES": {},
            "RESIDUES": self.residues,
            "MODS": self.mods,
            "SEPARATORS": self.separators,
        }
        n_rules = self.raw_rules["LMSD_CLASSES"]
        for c in self.supported_lmsd_classes:
            c_info = n_rules[c]
            seg_name_lst = c_info.get("ORDER", [])
            for seg_name in seg_name_lst:
                if seg_name in c_info:
                    seg_pattern_lst = c_info[seg_name]["INPUT"]
                    seg_pat_dct = {}
                    for seg_pattern in seg_pattern_lst:
                        seg_pat_dct[re.compile(seg_pattern)] = c_info[seg_name][
                            "OUTPUT"
                        ]
                    c_info[seg_name] = seg_pat_dct
            rules["LMSD_CLASSES"][c] = c_info
        self.rules = rules
        return rules
