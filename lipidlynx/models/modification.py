# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import re
from typing import Dict, List, Union

from jsonschema import Draft7Validator

from lipidlynx.controllers.general_functions import check_json, get_abs_path

from lipidlynx.models.defaults import (
    cv_elements_info,
    elem_nominal_info,
    lynx_schema,
    mod_level_lst,
)
from lipidlynx.models.log import logger
from lipidlynx.models.patterns import (
    mod_rgx,
    mod_db_rgx,
    mod_lv0_delta_rgx,
    mod_lv1_elem_rgx,
    mod_lv2_types_rgx,
    mod_lv3_position_rgx,
    mod_lv4_position_rgx,
)


class Modifications(object):
    def __init__(self, mod_code: str, db: int = 0):

        self.mod_code = mod_code.strip("<>")
        self.schema = "lynx_mod"
        with open(get_abs_path(lynx_schema[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(json.load(s_obj))

        self.db_count = db

        self.mod_info = self.__post_init__()
        self.mod_level = self.__identify_level__()
        self.sum_mod_info = self.get_sum_info()

        self.mod_id = self.sum_mod_info.get("mod_id", "")
        self.mod_linked_ids = self.sum_mod_info.get("mod_linked_ids", {})
        self.mod_list = self.sum_mod_info.get("mod_list", {})

        logger.info(
            f"Level {self.mod_level:3s} modification created from: {self.mod_code}"
        )

    def __parse_db__(self, db_str: str) -> dict:

        if mod_db_rgx.match(db_str):
            db_position_lst = db_str.strip("{}").split(",")
            if len(db_position_lst) != self.db_count:
                self.db_count = len(db_position_lst)

            db_dct = {"cv": "DB", "count": self.db_count}
            db_positions_info = self.get_positions(db_position_lst)
            if db_positions_info:
                db_dct.update(db_positions_info)
            return db_dct
        else:
            raise ValueError(f"Cannot parse DB info from string: {db_str}")

    def __parse_mod__(self, mod_str: str) -> dict:
        if mod_lv0_delta_rgx.match(self.mod_code):
            lv0_seg_lst = mod_lv0_delta_rgx.match(self.mod_code).groups()
            lv0_seg_lst = [s.strip(",") for s in lv0_seg_lst]
            delta_i = 0
            for s in lv0_seg_lst:
                delta_i += int(s)
            mod_dct = {
                "cv": "Delta",
                "count": delta_i,
                "positions": [],
                "positions_type": [],
            }
            return mod_dct
        else:
            mod_match = mod_lv3_position_rgx.match(mod_str)
            if mod_match:
                mod_match_dct = mod_match.groupdict()
                cv = mod_match_dct["cv"]
                count = int(
                    mod_match_dct.get("count") or 1
                )  # count can be '', set to 1 by default
                positions = mod_match_dct.get("positions", "")

                mod_dct = {"cv": cv, "count": count}

                if positions:
                    positions_info = self.get_positions(positions)
                    mod_dct.update(positions_info)
                return mod_dct
            else:
                if "+" in mod_str or "-" in mod_str:
                    mod_match = mod_lv1_elem_rgx.match(mod_str)
                else:
                    mod_match = mod_lv2_types_rgx.match(mod_str)
                if mod_match:
                    mod_match_dct = mod_match.groupdict()
                    cv = mod_match_dct["cv"]
                    count = int(mod_match_dct.get("count") or 1)
                    mod_dct = {
                        "cv": cv,
                        "count": count,
                        "positions": [],
                        "positions_type": [],
                    }
                    return mod_dct
                else:
                    raise ValueError(
                        f"Cannot parse modification info from string: {mod_str}"
                    )

    def __identify_level__(self) -> str:
        max_level = -1.0
        mod_level = -1
        if mod_lv0_delta_rgx.match(self.mod_code):
            mod_level = 0
        if mod_lv1_elem_rgx.match(self.mod_code):
            mod_level = 1
        if mod_lv2_types_rgx.match(self.mod_code):
            mod_level = 2
        for mod_seg in self.mod_info:
            mod_position_info = mod_seg.get("positions_info", None)
            if mod_seg.get("positions", None):
                if mod_level <= 3:
                    mod_level = 3
            if mod_position_info:
                for p_info_str in mod_position_info:
                    if re.match(r",?\d{1,2}[RS]", p_info_str):
                        if mod_level < 4:
                            mod_level = 4
        if self.mod_info:
            if self.mod_info[0]["cv"] == "DB":
                db_positions = self.mod_info[0].get("positions", None)
                db_positions_info = self.mod_info[0].get("positions_info", None)
                if db_positions:
                    max_level = max(max_level, mod_level + 0.1)
                if db_positions_info:
                    for db_info in db_positions_info:
                        if "E" in db_info or "Z" in db_info:
                            max_level = max(max_level, mod_level + 0.2)
                            break

        mod_level = f"{max(max_level, mod_level):.1f}"
        if mod_level.endswith(".0"):
            mod_level = mod_level[0]

        if float(mod_level) >= 0:
            return mod_level
        else:
            raise ValueError(
                f"Cannot parse the level of modification from string: {self.mod_code}"
            )

    def __post_init__(self):
        mod_info = []
        mod_code_match = mod_rgx.findall(self.mod_code)
        if mod_code_match:
            mod_lst = [
                m_seg[0].strip(",")
                for m_seg in mod_code_match
                if m_seg and m_seg[0].strip(",")
            ]
            if mod_lst:
                if mod_lst[0].startswith("{"):
                    db_dct = self.__parse_db__(mod_lst[0])
                    if db_dct:
                        mod_info.append(db_dct)
                    mod_lst = mod_lst[1:]

                for mod in mod_lst:
                    mod_dct = self.__parse_mod__(mod)
                    if mod_dct:
                        mod_info.append(mod_dct)

        return mod_info

    @staticmethod
    def __add_angle_brackets__(mod_str: str) -> str:
        return f"<{mod_str}>"

    @staticmethod
    def get_positions(positions: Union[str, List[str]]) -> dict:
        position_info = {"positions": [], "positions_info": []}
        if isinstance(positions, str):
            positions = positions.split(",")
        for position_seg in positions:
            position_match = mod_lv4_position_rgx.match(position_seg)
            if position_match:
                position_match_dct = position_match.groupdict()
                position = int(position_match_dct["position"])
                # position_type = position_match_dct.get('position_type', '')
                position_info["positions"].append(position)
                position_info["positions_info"].append(position_seg.strip(","))

        return position_info

    def get_elem_info(self):
        sum_mod_elem_dct = {}
        for mod in self.mod_info:
            mod_cv = mod["cv"]
            mod_count = mod["count"]
            mod_elem_dct = cv_elements_info[mod_cv]
            if mod_cv != "DB":
                for elem in mod_elem_dct:
                    if elem in sum_mod_elem_dct:
                        sum_mod_elem_dct[elem] += mod_count * mod_elem_dct[elem]
                    else:
                        sum_mod_elem_dct[elem] = mod_count * mod_elem_dct[elem]

        return sum_mod_elem_dct

    def get_sum_info(self):
        linked_ids_lst = self.to_all_levels()
        mod_id = linked_ids_lst.get(self.mod_level, "")
        if mod_id:
            if mod_id == self.mod_code:
                sum_mod_info_dct = {
                    "_version": lynx_schema.get("_version", "0.1"),
                    "mod_id": self.mod_code,
                    "mod_level": self.mod_level,
                    "mod_linked_ids": self.to_all_levels(),
                    "mod_list": self.mod_info,
                }
            else:
                sum_mod_info_dct = {
                    "_version": lynx_schema.get("_version", "0.1"),
                    "mod_id": mod_id,
                    "input_mod_id": self.mod_code,
                    "mod_level": self.mod_level,
                    "mod_linked_ids": self.to_all_levels(),
                    "mod_list": self.mod_info,
                }
        else:
            raise ValueError(
                f"Cannot format modification code to level {self.mod_level} "
                f"from input: {self.mod_code}"
            )

        return sum_mod_info_dct

    def to_json(self):
        mod_json_str = json.dumps(self.sum_mod_info)
        if check_json(validator=self.validator, json_obj=json.loads(mod_json_str)):
            return mod_json_str
        else:
            raise Exception(f"Schema test FAILED. Schema {self.schema}")

    def to_mass_shift(self, angle_brackets: bool = True) -> str:
        if float(self.mod_level) > 0:
            sum_mod_elem_dct = self.get_elem_info()
            delta = 0
            for elem in sum_mod_elem_dct:
                if elem in elem_nominal_info:
                    delta += elem_nominal_info[elem] * sum_mod_elem_dct[elem]

            mod_str = f"{delta:+}"
            if angle_brackets:
                mod_str = self.__add_angle_brackets__(mod_str)
            return mod_str
        elif float(self.mod_level) == 0:
            lv0_match = mod_lv0_delta_rgx.match(self.mod_code)
            lv0_groups = []
            if lv0_match:
                lv0_matched_groups = lv0_match.groups()
                lv0_groups = [int(i.strip(",")) for i in lv0_matched_groups]
            if lv0_groups:
                mod_str = f"{sum(lv0_groups):+}"
                if angle_brackets:
                    mod_str = self.__add_angle_brackets__(mod_str)
            else:
                mod_str = ""
            return mod_str
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 0, received: {self.mod_level}"
            )

    def to_elem_shift(self, angle_brackets: bool = True) -> str:
        if float(self.mod_level) >= 1:
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

            mod_str = f'{",".join(mod_str_lst)}'
            if angle_brackets:
                mod_str = self.__add_angle_brackets__(mod_str)
            return mod_str
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 1, received: {self.mod_level}"
            )

    def to_mod_count(self, angle_brackets: bool = True) -> str:
        if float(self.mod_level) >= 2:
            mod_output_lst = []
            for mod in self.mod_info:
                if mod["cv"] != "DB":
                    if mod["count"] > 1:
                        mod_output_lst.append(f'{mod["count"]}{mod["cv"]}')
                    else:
                        mod_output_lst.append(f'{mod["cv"]}')
            mod_str = f'{",".join(mod_output_lst)}'
            if angle_brackets:
                mod_str = self.__add_angle_brackets__(mod_str)
            return mod_str
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 2, received: {self.mod_level}"
            )

    def to_mod_position(
        self,
        angle_brackets: bool = True,
        db_position: bool = False,
        db_position_info=False,
    ) -> str:
        if db_position_info and not db_position:
            db_position = True
        if float(self.mod_level) >= 3:
            mod_output_lst = []
            for mod in self.mod_info:
                positions_lst = [str(i) for i in mod["positions"]]
                positions = "{" + ",".join(positions_lst) + "}"
                if mod["cv"] == "DB":
                    if db_position and db_position_info:  # mod level 3.2
                        mod_output_lst.append(
                            "{" + ",".join(mod["positions_info"]) + "}"
                        )
                    if db_position and not db_position_info:  # mod level 3.1
                        mod_output_lst.append(positions)
                else:
                    if mod["count"] > 1:
                        mod_output_lst.append(f'{mod["count"]}{mod["cv"]}{positions}')
                    else:
                        mod_output_lst.append(f'{mod["cv"]}{positions}')
            mod_str = f'{",".join(mod_output_lst)}'
            if angle_brackets:
                mod_str = self.__add_angle_brackets__(mod_str)
            return mod_str
        else:
            raise ValueError(
                f"Modification level not fit. required level >= 3, received: {self.mod_level}"
            )

    def to_mod_position_type(
        self,
        angle_brackets: bool = True,
        db_position: bool = False,
        db_position_info: bool = False,
    ) -> str:
        if db_position_info and not db_position:
            db_position = True

        if float(self.mod_level) >= 4:
            mod_output_lst = []
            for mod in self.mod_info:
                positions = "{" + ",".join(mod["positions_info"]) + "}"
                if mod["cv"] == "DB":
                    if db_position and db_position_info:  # mod level 4.2
                        mod_output_lst.append(positions)
                    if db_position and not db_position_info:  # mod level 4.1
                        positions_lst = [str(i) for i in mod["positions"]]
                        mod_output_lst.append("{" + ",".join(positions_lst) + "}")
                else:
                    if mod["count"] > 1:
                        mod_output_lst.append(f'{mod["count"]}{mod["cv"]}{positions}')
                    else:
                        mod_output_lst.append(f'{mod["cv"]}{positions}')
            mod_str = f'{",".join(mod_output_lst)}'
            if angle_brackets:
                mod_str = self.__add_angle_brackets__(mod_str)
            return mod_str
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
        if level in mod_level_lst:
            if level == "0":
                mod_str = self.to_mass_shift(angle_brackets=angle_brackets)
            elif level == "1":
                mod_str = self.to_elem_shift(angle_brackets=angle_brackets)
            elif level == "2":
                mod_str = self.to_mod_count(angle_brackets=angle_brackets)
            elif level == "3":
                mod_str = self.to_mod_position(
                    angle_brackets=angle_brackets, db_position=False
                )
            elif level == "3.1":
                mod_str = self.to_mod_position(
                    angle_brackets=angle_brackets, db_position=True
                )
            elif level == "3.2":
                mod_str = self.to_mod_position(
                    angle_brackets=angle_brackets,
                    db_position=True,
                    db_position_info=True,
                )
            elif level == "4":
                mod_str = self.to_mod_position_type(
                    angle_brackets=angle_brackets,
                    db_position=False,
                    db_position_info=False,
                )
            elif level == "4.1":
                mod_str = self.to_mod_position_type(
                    angle_brackets=angle_brackets,
                    db_position=True,
                    db_position_info=False,
                )
            elif level == "4.2":
                mod_str = self.to_mod_position_type(
                    angle_brackets=angle_brackets,
                    db_position=True,
                    db_position_info=True,
                )
            else:
                raise ValueError(f"Currently not supported modification level: {level}")

        return mod_str

    def to_all_levels(
        self, angle_brackets: bool = True, as_list: bool = False
    ) -> Union[Dict[str, str], List[str]]:
        all_levels_dct = {}
        all_levels_lst = []
        if self.mod_level in mod_level_lst:
            mod_idx = mod_level_lst.index(self.mod_level)
            output_levels_lst = mod_level_lst[: mod_idx + 1]
        else:
            raise ValueError(
                f"Currently not supported modification level: {self.mod_level}"
            )
        if self.mod_level == "4":
            output_levels_lst.remove("3.2")
            output_levels_lst.remove("3.1")
        elif self.mod_level == "4.1":
            output_levels_lst.remove("3.2")
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
        r"<-18>",
        r"<+46>",
        r"<+3O,-2H>",
        r"<2OH,Ke>",
        r"<2OH{8,11},Ke{14}>",
        r"<{5,9,12,15},2OH{8,11},Ke{14}>",
        r"<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>",
        r"<2OH{8R,11S},Ke{14}>",
        r"<{5,9,12,15},2OH{8R,11S},Ke{14}>",
        r"<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>",
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
