# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os
from typing import Dict, List, Union

from jsonschema import Draft7Validator, RefResolver


from lynx.models.cv import CV
from lynx.models.defaults import (
    api_version,
    cv_elements_info,
    elem_nominal_info,
    lynx_schema_cfg,
    core_schema,
    core_schema_path,
    mod_db_level_lst,
    default_cv_file
)
from lynx.utils.file_readers import get_abs_path
from lynx.utils.log import logger
from lynx.utils.toolbox import check_json


class Mods(object):

    def __init__(self, mod_info: dict, db: int = 0, schema="lynx_mod", cv_file: str = default_cv_file):
        if os.path.isfile(cv_file):
            pass
        else:
            cv_file = default_cv_file
        self.cv = CV(cv_file).info
        self.mod_info = mod_info
        self.schema = "lynx_mod"
        self.type = "Modification"
        self.mod_level = self.mod_info.get("MOD_LEVEL", 0)
        with open(get_abs_path(lynx_schema_cfg[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(
                json.load(s_obj),
                resolver=RefResolver(
                    f"file://{core_schema_path}", referrer=core_schema
                ),
            )

        self.db_count = db

        self.mod_info = self.__post_init__()

        self.sum_mod_info = self.get_sum_info()

        self.mod_id = self.sum_mod_info.get("id", "")
        self.mod_linked_ids = self.sum_mod_info.get("linked_ids", {})
        self.mod_list = self.sum_mod_info.get("info", {})

        # logger.info(
        #     f"Level {self.mod_level:3s} modification created from: {self.mod_code}"
        # )

    @staticmethod
    def __add_angle_brackets__(mod_str: str) -> str:
        if mod_str and mod_str not in ["+0", "-0", "0"]:
            return f"<{mod_str}>"
        else:
            return ""

    def get_sum_elements_shift(self):
        sum_mod_elem_dct = {}
        for mod in self.mod_info:
            mod_cv = mod["MOD_CV"]
            mod_count = mod["MOD_COUNT"]
            mod_elem_dct = mod["MOD_ELEMENTS"]
            if mod_cv != "DB":
                for elem in mod_elem_dct:
                    if elem in sum_mod_elem_dct:
                        sum_mod_elem_dct[elem] += mod_count * mod_elem_dct[elem]
                    else:
                        sum_mod_elem_dct[elem] = mod_count * mod_elem_dct[elem]

        return sum_mod_elem_dct

    def get_sum_info(self):
        linked_ids = self.to_all_levels()
        mod_id = linked_ids.get(self.mod_level, "")
        if mod_id:
            if mod_id == self.mod_code:
                sum_mod_info_dct = {
                    "api_version": api_version,
                    "type": self.type,
                    "id": self.mod_code,
                    "level": self.mod_level,
                    "linked_ids": self.to_all_levels(),
                    "info": self.mod_info,
                }
            else:
                sum_mod_info_dct = {
                    "api_version": api_version,
                    "type": self.type,
                    "id": mod_id,
                    # "input_name": self.mod_code,
                    "level": self.mod_level,
                    "linked_ids": self.to_all_levels(),
                    "info": self.mod_info,
                }
        else:
            raise ValueError(
                f"Cannot format_mods modification code to level {self.mod_level} "
                f"from input: {self.mod_code}"
            )

        return sum_mod_info_dct

    def to_json(self):
        mod_json_str = json.dumps(self.sum_mod_info)

        if check_json(validator=self.validator, json_obj=json.loads(mod_json_str)):
            return mod_json_str
        else:
            raise Exception(f"Schema test FAILED. Schema {self.schema}")

    def to_mass_shift(self) -> str:
        if float(self.mod_level) > 1:
            sum_mod_elem_dct = self.get_elem_info()
            delta = 0
            for elem in sum_mod_elem_dct:
                if elem in elem_nominal_info:
                    delta += elem_nominal_info[elem] * sum_mod_elem_dct[elem]
            return f"{delta:+}"
        elif float(self.mod_level) == 1:
            lv0_match = mod_lv0_delta_rgx.match(self.mod_code)
            lv0_groups = []
            if lv0_match:
                lv0_matched_groups = lv0_match.groups()
                lv0_groups = [int(i.strip(",")) for i in lv0_matched_groups]
            if lv0_groups:
                mod_str = f"{sum(lv0_groups):+}"
            else:
                mod_str = ""
            return mod_str
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 1, received: {self.mod_level}"
            )

    def to_elem_shift(self) -> str:
        if float(self.mod_level) >= 2:
            sum_mod_elem_dct = self.get_elem_info()
            mod_str_lst = []
            mod_elem_lst = ["O", "N", "S", "H", "Na"]
            for mod_elem in mod_elem_lst:
                if mod_elem in sum_mod_elem_dct:
                    mod_elem_count = f"{sum_mod_elem_dct[mod_elem]:+}"
                    if mod_elem_count == "+1":
                        mod_elem_count = "+"
                    elif mod_elem_count == "-1":
                        mod_elem_count = "-"
                    mod_str_lst.append(f"{mod_elem_count}{mod_elem}")

            return f'{",".join(mod_str_lst)}'
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 2, received: {self.mod_level}"
            )

    def to_mod_count(self) -> str:
        if float(self.mod_level) >= 3:
            mod_output_lst = []
            for mod in self.mod_info:
                if mod["cv"] != "DB":
                    if mod["count"] > 1:
                        mod_output_lst.append(f'{mod["count"]}{mod["cv"]}')
                    else:
                        mod_output_lst.append(f'{mod["cv"]}')
            return f'{",".join(mod_output_lst)}'
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 3, received: {self.mod_level}"
            )

    def to_mod_position(self) -> str:
        if float(self.mod_level) >= 4:
            mod_output_lst = []
            for mod in self.mod_info:
                positions_lst = [str(i) for i in mod["positions"]]
                positions = "{" + ",".join(positions_lst) + "}"
                if mod["cv"] == "DB":
                    pass
                else:
                    if mod["count"] > 1:
                        mod_output_lst.append(f'{mod["count"]}{mod["cv"]}{positions}')
                    else:
                        mod_output_lst.append(f'{mod["cv"]}{positions}')
            return f'{",".join(mod_output_lst)}'
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 4, received: {self.mod_level}"
            )

    def to_mod_position_type(self) -> str:

        if float(self.mod_level) >= 5:
            mod_output_lst = []
            for mod in self.mod_info:
                positions = "{" + ",".join(mod["positions_info"]) + "}"
                if mod["cv"] == "DB":
                    pass
                else:
                    if mod["count"] > 1:
                        mod_output_lst.append(f'{mod["count"]}{mod["cv"]}{positions}')
                    else:
                        mod_output_lst.append(f'{mod["cv"]}{positions}')

            return f'{",".join(mod_output_lst)}'
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 4, received: {self.mod_level}"
            )

    def to_mod_level(
        self, level: Union[int, float, str] = 0, angle_brackets: bool = True
    ) -> str:
        mod_str = ""
        if not isinstance(level, str):
            level = str(level)
        if float(level) > float(self.mod_level):
            raise ValueError(
                f'Cannot convert to higher level than the mod_level "{self.mod_level}". Input:{level}'
            )
        if level in mod_db_level_lst and not level.startswith("0"):
            if level.startswith("1"):
                mod_str = self.to_mass_shift()
            elif level.startswith("2"):
                mod_str = self.to_elem_shift()
            elif level.startswith("3"):
                mod_str = self.to_mod_count()
            elif level.startswith("4"):
                mod_str = self.to_mod_position()
            elif level.startswith("5"):
                mod_str = self.to_mod_position_type()
            else:
                raise ValueError(f"Currently not supported modification level: {level}")
        if self.mod_info:
            for mod in self.mod_info:
                if mod["cv"] == "DB" and mod.get("positions", None):
                    if level.endswith(".1"):
                        positions_lst = [str(i) for i in mod["positions"]]
                        positions = "{" + ",".join(positions_lst) + "}"
                        mod_str = ",".join([positions, mod_str])
                    if level.endswith(".2"):
                        mod_str = ",".join(
                            ["{" + ",".join(mod["positions_info"]) + "}", mod_str]
                        )
        mod_str = mod_str.strip(",")
        if mod_str and angle_brackets:
            mod_str = self.__add_angle_brackets__(mod_str)

        return mod_str

    def to_all_levels(
        self, angle_brackets: bool = True, as_list: bool = False
    ) -> Union[Dict[str, str], List[str]]:
        all_levels_dct = {}
        all_levels_lst = []
        if self.mod_level in mod_db_level_lst:
            mod_idx = mod_db_level_lst.index(self.mod_level)
            output_levels_lst = mod_db_level_lst[: mod_idx + 1]
        else:
            raise ValueError(f"Modification level not supported: {self.mod_level}")
        if len(self.mod_level) == 1:
            output_levels_lst = [
                out_lv for out_lv in output_levels_lst if len(out_lv) == 1
            ]
        elif self.mod_level.endswith(".1"):
            output_levels_lst = [
                out_lv
                for out_lv in output_levels_lst
                if len(out_lv) == 1 or out_lv.endswith(".1")
            ]
        else:
            pass
        if as_list:
            for level in output_levels_lst:
                all_levels_lst.append(
                    self.to_mod_level(level, angle_brackets=angle_brackets)
                )
            all_levels_info = all_levels_lst
        else:
            for level in output_levels_lst:
                all_levels_dct[level] = self.to_mod_level(
                    level, angle_brackets=angle_brackets
                )
            all_levels_info = all_levels_dct

        return all_levels_info


if __name__ == "__main__":

    mod_code_lst = [
        # r"<-18>",
        # r"<+46>",
        # r"<+3O,-2H>",
        # r"<2OH,Ke>",
        # r"<2OH{8,11},Ke{14}>",
        r"<{5,9,12,15},2OH{8,11},oxo{14}>",
        r"<{5Z,9E,12E,15E},2OH{8,11},oxo{14}>",
        r"<2OH{8R,11S},oxo{14}>",
        r"<{5,9,12,15},2OH{8R,11S},oxo{14}>",
        r"<{5Z,9E,12E,15E},2OH{8R,11S},oxo{14}>",
        r"<{9}>",
        r"<{9,12}>",
        r"<{9Z}>",
        r"<{9Z,11E}>",
    ]

    for usr_mod_code in mod_code_lst:
        logger.info(f"Test mod str: {usr_mod_code}")

        mod_obj = Modifications(usr_mod_code)

        # print('to mass shift', mod_obj.to_mass_shift())
        # print('to elem shift', mod_obj.to_elem_shift())
        # print('to mod_count', mod_obj.to_mod_count())
        # print('to mod_position', mod_obj.to_mod_position())
        # print('to mod_position_type', mod_obj.to_mod_position_type())
        # print('to mass shift', mod_obj.to_mass_shift(angle_brackets=False))
        # print('to elem shift', mod_obj.to_elem_shift(angle_brackets=False))
        # print('to mod_count', mod_obj.to_mod_count(angle_brackets=False))
        # print('to mod_position', mod_obj.to_mod_position(angle_brackets=False))
        # print('to mod_position_type', mod_obj.to_mod_position_type(angle_brackets=False))
        logger.debug(f"to all levels: {mod_obj.to_all_levels()}")
        logger.debug(
            f"to all levels without <> : {mod_obj.to_all_levels(angle_brackets=False)}"
        )

        mod_json = mod_obj.to_json()
        logger.debug(mod_json)

    logger.info("FINISHED")
