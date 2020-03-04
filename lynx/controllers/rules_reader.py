import json
import os
import re
from typing import Union

from lynx.models.log import logger
from lynx.controllers.general_functions import get_abs_path, load_folder


def js_reader(file):
    file = get_abs_path(file)
    with open(file) as file_obj:
        js_obj = json.load(file_obj)
        print(js_obj, type(js_obj))
        return js_obj


class Rules(object):
    """
    Read rules from json file and build corresponding regular expressions
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
        self._date = self.raw_rules.get("_DATE", 20200214)
        self._authors = self.raw_rules.get("_AUTHORS", ["example@uni-example.de"])
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
        logger.info(
            f"Load rule settings: for {self.source}\n"
            f"Last modified: {self._date}\n"
            f"Authors: {self._authors}"
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

    def __build__(self, lipid_classes: list, group: str):
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

                rules[c] = {"GROUPS": seg_lst, "PATTERN": pattern_str}
        return rules

    def build(self) -> dict:
        sum_rules = self.__build__(self.supported_residues, "RESIDUES")
        self.rules = sum_rules
        sum_rules.update(self.__build__(self.supported_classes, "LIPID_CLASSES"))
        self.rules = sum_rules
        return sum_rules

    def validate(self) -> bool:
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


def build_rules(folder: str) -> dict:

    sum_rules = {}
    file_path_lst = load_folder(folder, file_type=".json")
    logger.debug(f"Fund JSON config files: \n {file_path_lst}")

    for f in file_path_lst:
        temp_rules = Rules(f)
        idx_lst = [os.path.basename(f)] + temp_rules.source
        sum_rules["#".join(idx_lst)] = temp_rules

    logger.info(sum_rules)

    return sum_rules


if __name__ == "__main__":
    file_js = r"../configurations/rules/input/MS-DIAL.json"
    usr_rules = Rules(file_js)
    out_dct = usr_rules.build()
    usr_rules.validate()

    js_folder = r"../configurations/rules/input"
    build_rules(js_folder)

    logger.info("Finished.")
