# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
from typing import Union

from lynx.models.log import logger
from lynx.controllers.general_functions import js_reader


class InputRules(object):
    """
    Read input rules from json file and build corresponding regular expressions
    """

    def __init__(self, rules: Union[str, dict]):
        if isinstance(rules, dict):
            pass
        elif isinstance(rules, str):
            rules = js_reader(rules)
        else:
            raise TypeError
        self.raw_rules = rules
        self.source = self.raw_rules["SOURCE"]
        self.date = self.raw_rules.get("_DATE", 20200214)
        self.authors = self.raw_rules.get("_AUTHORS", ["example@uni-example.de"])
        self.supported_residues = list(self.raw_rules["RESIDUES"].keys())
        self.supported_classes = list(self.raw_rules["LIPID_CLASSES"].keys())
        self.separators = self.raw_rules["SEPARATORS"]
        self.rules = {}
        self.is_structure_valid = self.__check__()
        self.__replace_fields__()
        if self.is_structure_valid:
            pass
        else:
            raise ValueError(
                "Cannot find valid sections for key [SOURCE, RESIDUES, LIPID_CLASSES]"
            )
        self.rules = self.build()
        self.is_validated = self.validate()
        logger.info(
            f"Load input rule: {self.source}\n"
            f"Last modified: {self.date}\n"
            f"Authors: {self.authors}"
        )

    def __replace_fields__(self):
        for sep in self.separators:
            for res in self.supported_residues:
                res_sep = self.raw_rules["RESIDUES"][res].get(sep, None)
                if res_sep == sep:
                    self.raw_rules["RESIDUES"][res][sep] = self.separators.get(
                        sep, "_|/"
                    )
            for lc in self.supported_classes:
                lc_sep = self.raw_rules["LIPID_CLASSES"][lc].get(sep, None)
                if lc_sep == sep:
                    self.raw_rules["LIPID_CLASSES"][lc][sep] = self.separators.get(
                        sep, r"\s"
                    )

    def __check__(self):
        is_structure_valid = False
        if (
            isinstance(self.source, list)
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
                pattern_str = ""
                seg_lst = []
                for seg in order_lst:
                    seg_str = temp_c_dct.get(seg, "")
                    if (
                        c in self.supported_classes
                        and seg_str in self.supported_residues
                    ):
                        seg_str = self.rules[seg_str].get("PATTERN", "")
                        if seg_str:
                            pattern_str += seg_str
                            seg_lst.append(seg)
                        else:
                            raise KeyError(f"Rule pattern {seg_str} not defined!")
                    else:
                        pattern_str += f"(?P<{seg}>{seg_str})"
                        seg_lst.append(seg)

                    if seg in optional_lst:
                        pattern_str += "?"

                rules[c] = {
                    "GROUPS": seg_lst,
                    "PATTERN": pattern_str,
                    "MATCH": re.compile(pattern_str),
                    "CLASS": temp_c_dct.get("CLASS", c),
                    "LMSD_CLASSES": temp_c_dct.get("LMSD_CLASSES", [c]),
                    "RESIDUES_SEPARATOR": self.separators.get(
                        "RESIDUES_SEPARATOR", "_|/"
                    ),
                }
        return rules

    def build(self) -> dict:
        sum_rules = self.__build__(self.supported_residues, "RESIDUES")
        self.rules = sum_rules
        sum_rules.update(self.__build__(self.supported_classes, "LIPID_CLASSES"))
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
        for c in self.supported_classes:
            temp_c_dct = self.raw_rules.get("LIPID_CLASSES", {}).get(c, {})
            pattern_dct = self.rules.get(c, "")
            if isinstance(pattern_dct, dict):
                pattern_str = pattern_dct.get("PATTERN", None)
            else:
                pattern_str = None
            if pattern_str and temp_c_dct:

                max_res_count = self.raw_rules["LIPID_CLASSES"][c].get(
                    "MAX_RESIDUES", 1
                )
                logger.debug(f"Validating {c} pattern: {pattern_str}")
                c_pattern = re.compile(pattern_str)
                test_lst = temp_c_dct.get("EXAMPLE", [])
                test_dct = {}
                for test_str in test_lst:
                    fit_class = False
                    class_rgx = re.compile(pattern_dct["CLASS"])
                    class_search = class_rgx.search(test_str)
                    if class_search:
                        logger.debug(f"{test_str} fit to class {c}")
                    else:
                        logger.error(f"{test_str} NOT fit to class {c}")
                    c_match = c_pattern.match(test_str)
                    if c_match:
                        c_matched_res_dct = c_match.groupdict()
                        logger.debug(f"Check: {test_str}")
                        logger.debug(c_matched_res_dct)
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
                                        logger.debug(f"Check: {test_str} # {res}")
                                        logger.debug(res_pat.match(res).groupdict())
                                    else:
                                        valid_lst.append(False)
                                        test_dct[f"{test_str} # {res}"] = "failed"
                                        logger.error(f"{test_str} # {res} failed")
                            else:
                                raise ValueError(
                                    f"Number of residues exceed the defined limit: "
                                    f"{max_res_count} for Lipid_Class: {c} from source: {self.source}"
                                )
                        else:
                            valid_lst.append(True)
                            test_dct[test_str] = "passed"
                    else:
                        valid_lst.append(False)
                        test_dct[test_str] = "failed"
                        logger.error(f"{test_str} failed")
                logger.debug(test_dct)
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
            rules = js_reader(rules)
        else:
            raise TypeError
        self.raw_rules = rules
        self.date = self.raw_rules.get("_DATE", 20200214)
        self.authors = self.raw_rules.get("_AUTHORS", ["example@uni-example.de"])
        self.nomenclature = self.raw_rules.get("NOMENCLATURE", "LipidLynxX")
        self.supported_lmsd_classes = list(self.raw_rules["LMSD_CLASSES"].keys())
        self.separators = self.raw_rules["SEPARATORS"]
        self.rules = self.build()
        self.is_structure_valid = self.__check__()
        logger.info(
            f"Load rules for nomenclature {self.nomenclature}\n"
            f"Last modified: {self.date}\n"
            f"Authors: {self.authors}"
        )

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
        rules = {"SEPARATORS": self.separators, "LMSD_CLASSES": {}}
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
