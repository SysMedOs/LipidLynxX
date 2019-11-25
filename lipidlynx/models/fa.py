import json
import pandas as pd
import re
from typing import Dict, List, Union

from jsonschema import Draft7Validator

from lipidlynx.controllers.Logger import logger

from lipidlynx.models.defaults import (
    lynx_schema,
    cv_elements_info,
    elem_nominal_info,
    lipid_level_lst,
)

from lipidlynx.models.modification import Modifications
from lipidlynx.models.patterns import fa_rgx, mod_rgx

from lipidlynx.controllers.GeneralFunctions import get_abs_path


class FattyAcid(object):
    def __init__(self, lipid_code: str, db: int = 0):

        self.lipid_code = lipid_code.strip("FA")

        with open(get_abs_path(lynx_schema["lynx_fa"]), "r") as s_obj:
            self.validator = Draft7Validator(json.load(s_obj))

        self.fa_info_dct = self.__post_init__()
        self.fa_info_dct["fa_id"] = self.lipid_code

        self.mod_info = self.fa_info_dct.get("mod_obj", None)
        if db and db < self.fa_info_dct["fa_segments"]["db"]:
            self.db_count = db
        else:
            self.db_count = self.fa_info_dct["fa_segments"]["db"]
        self.fa_level = self.fa_info_dct.get("fa_level", "")
        self.fa_linked_ids = self.fa_info_dct.get("fa_linked_ids", {})
        logger.info(
            f"Level {self.fa_level:4s} FattyAcid created from: {self.lipid_code}"
        )

    def __post_init__(self):

        fa_info_dct = {}
        fa_match = fa_rgx.match(self.lipid_code)

        if fa_match:
            fa_matched_dct = fa_match.groupdict()
            fa_info_dct["fa_segments"] = fa_matched_dct
            mod_code = fa_matched_dct.get("mod", None)

            fa_seg_str = f'{fa_matched_dct.get("link", "")}{fa_matched_dct.get("c", "")}:{fa_matched_dct.get("db", "")}'
            if fa_seg_str.lower().startswith('none'):
                fa_seg_str = fa_seg_str[4:]
            if not fa_seg_str[:2] in ['O-', 'P-', 'FA']:
                fa_seg_str = f'FA{fa_seg_str}'
            if mod_code:
                mod_obj = Modifications(mod_code)
                fa_info_dct["mod_obj"] = mod_obj
                fa_info_dct["mod_info"] = mod_obj.sum_mod_info
                fa_info_dct["mod_id"] = mod_obj.mod_id
                fa_info_dct["fa_level"] = f"S{mod_obj.mod_level}"
                fa_linked_ids_dct = {}
                mod_linked_ids_dct = mod_obj.mod_linked_ids  # type: dict
                for lv in mod_linked_ids_dct:
                    fa_linked_ids_dct[
                        f"S{lv}"
                    ] = f"{fa_seg_str}{mod_linked_ids_dct[lv]}"
                fa_info_dct["fa_linked_ids"] = fa_linked_ids_dct
                logger.debug(f"{fa_info_dct}")
            return fa_info_dct
        else:
            raise ValueError(f"Cannot parse FA sting: {self.lipid_code}")

    def to_json(self):
        fa_lite_info_dct = self.fa_info_dct
        fa_lite_info_dct.pop("mod_obj", None)
        fa_json_str = json.dumps(fa_lite_info_dct)
        j_json = json.loads(fa_json_str)
        if self.validator.is_valid(j_json):
            logger.debug(f"Schema test PASSED. FA JSON created.")
            return fa_json_str
        else:
            errors = sorted(self.validator.iter_errors(j_json), key=lambda e: e.path)
            for e in errors:
                logger.error(e)
            raise Exception(f"Schema test FAILED.")


if __name__ == "__main__":

    mod_code_lst = [
        r"FA20:4<-18>",
        r"FA20:4<+46>",
        r"FA20:4<+3O,-2H>",
        r"FA20:4<2OH,Ke>",
        r"FA20:4<2OH{8,11},Ke{14}>",
        r"FA20:4<{5,9,12,15},2OH{8,11},Ke{14}>",
        r"FA20:4<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>",
        r"FA20:4<2OH{8R,11S},Ke{14}>",
        r"FA20:4<{5,9,12,15},2OH{8R,11S},Ke{14}>",
        r"FA20:4<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>",
    ]

    for usr_mod_code in mod_code_lst:
        logger.info(f"Test FA str: {usr_mod_code}")

        fa_obj = FattyAcid(usr_mod_code)

        # logger.debug(f"to all levels: {fa_obj.to_all_levels()}")
        # logger.debug(
        #     f"to all levels without <> : {fa_obj.to_all_levels(angle_brackets=False)}"
        # )
        #
        fa_json = fa_obj.to_json()
        logger.debug(fa_json)

    logger.info("FINISHED")
