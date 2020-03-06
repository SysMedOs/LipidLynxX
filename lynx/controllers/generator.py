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
    input_rules,
    output_rules,
)
from lynx.controllers.general_functions import seg_to_str
from lynx.controllers.params_loader import build_input_rules, build_output_rules
from lynx.controllers.parser import parse


def generate(lipid_name: str, rules: dict, rule: str) -> str:
    output_name = ""
    out_r_info = rules.get(rule, None)
    if out_r_info:
        parsed_info = parse(lipid_name, rules=input_rules)
        for p in parsed_info:
            p_info = parsed_info[p]
            logger.info(p_info)
            for in_r in p_info:
                r_info = p_info[in_r]  # type: dict
                lmsd_classes = r_info.get("LMSD_CLASSES", None)
                segments = r_info["SEGMENTS"]
                c_str = segments.get("CLASS", "")
                for c in lmsd_classes:
                    if c in out_r_info:
                        lmsd_patterns_lst = out_r_info[c].get("INPUT_PATTERNS")
                        for c_rgx in lmsd_patterns_lst:
                            logger.debug(f"Test {c_str} on LMSD: {c} using {c_rgx}")
                            c_searched = c_rgx.search(c_str)
                            if c_searched:
                                logger.info(f"{lipid_name} is LMSD: {c} by {in_r}")
    else:
        raise ValueError(f"Cannot find output rule: {rule}")

    logger.debug(parsed_info)
    return output_name


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

    examples = ["PE(P-16:0_18:2)"]

    for e in examples:
        o = generate(e, rules=output_rules, rule="LipidLynxX@20200214")
        logger.debug(e)
        logger.debug(o)

    logger.info("fin")
