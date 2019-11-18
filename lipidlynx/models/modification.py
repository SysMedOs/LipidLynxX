import json
import pandas as pd
import re
from typing import Union, List

from jsonschema import Draft7Validator

from lipidlynx.models.defaults import lynx_schema, cv_elements_info, logger, elem_nominal_info
from lipidlynx.models.patterns import mod_delta_rgx, mod_elem_rgx, mod_rgx, mod_db_rgx, mod_no_position_rgx, \
    mod_w_position_rgx, mod_position_rgx
from lipidlynx.controllers.GeneralFunctions import get_abs_path


class Modifications(object):

    def __init__(self, mod_code: str, db: int = 0):

        self.mod_code = mod_code.strip('<>')

        with open(get_abs_path(lynx_schema['lynx_mod']), 'r') as s_obj:
            self.validator = Draft7Validator(json.load(s_obj))

        self.db_count = db

        self.mod_level = self.__identify_level__()  # type: int
        print(self.mod_level)

        self.mod_info = self.__post_init__()

    def __parse_db__(self, db_str: str) -> dict:
        db_match = mod_db_rgx.match(db_str)
        if db_match:
            db_position_lst = db_match.groups()

            if len(db_position_lst) != self.db_count:
                self.db_count = len(db_position_lst)

            db_dct = {"cv": "DB", "count": self.db_count}
            db_positions_info = self.get_positions(db_position_lst)
            if db_positions_info:
                db_dct.update(db_positions_info)
            return db_dct
        else:
            raise ValueError(f'Cannot parse DB info from string: {db_str}')

    def __parse_mod__(self, mod_str: str) -> dict:
        if self.mod_level >= 1:
            mod_match = mod_w_position_rgx.match(mod_str)
            if mod_match:
                mod_match_dct = mod_match.groupdict()
                cv = mod_match_dct['cv']
                count = int(mod_match_dct.get('count') or 1)  # count can be '', set to 1 by default
                positions = mod_match_dct.get('positions', '')

                mod_dct = {"cv": cv, "count": count}

                if positions:
                    positions_info = self.get_positions(positions)
                    mod_dct.update(positions_info)
                return mod_dct
            else:
                mod_match = mod_no_position_rgx.match(mod_str)
                if mod_match:
                    mod_match_dct = mod_match.groupdict()
                    cv = mod_match_dct['cv']
                    count = int(mod_match_dct.get('count') or 1)
                    mod_dct = {"cv": cv, "count": count, 'positions': [], 'positions_type': []}
                    return mod_dct
                else:
                    raise ValueError(f'Cannot parse modification info from string: {mod_str}')
        else:
            mod_dct = {"cv": "Delta", "count": 1, 'positions': [], 'positions_type': []}
            return mod_dct

    def __identify_level__(self) -> int:
        mod_level = -1
        if mod_delta_rgx.match(self.mod_code):
            mod_level = 0
        if mod_elem_rgx.match(self.mod_code):
            print(mod_elem_rgx.match(self.mod_code).groups())
            mod_level = 1
        if mod_no_position_rgx.match(self.mod_code):
            mod_level = 2
        if mod_w_position_rgx.match(self.mod_code):
            mod_level = 3
        if mod_level >= 0:
            return mod_level
        else:
            raise ValueError(f'Cannot parse the level of modification from string: {self.mod_code}')

    def __post_init__(self):
        mod_info = []
        mod_code_match = mod_rgx.findall(self.mod_code)
        if mod_code_match:
            mod_lst = [m_seg[0].strip(',') for m_seg in mod_code_match if m_seg and m_seg[0].strip(',')]
            if mod_lst:
                if mod_lst[0].startswith('{'):
                    db_dct = self.__parse_db__(mod_lst[0])
                    if db_dct:
                        mod_info.append(db_dct)
                    mod_lst = mod_lst[1:]

                for mod in mod_lst:
                    print(mod)
                    mod_dct = self.__parse_mod__(mod)
                    if mod_dct:
                        mod_info.append(mod_dct)

        logger.info(f'Modification object created from string: {self.mod_code}')
        return mod_info

    @staticmethod
    def __add_angle_brackets__(mod_str: str) -> str:
        return f'<{mod_str}>'

    @staticmethod
    def get_positions(positions: Union[str, List[str]]) -> dict:
        position_info = {
            'positions': [],
            'positions_type': []
        }
        if isinstance(positions, str):
            positions = positions.split(',')
        for position_seg in positions:
            position_match = mod_position_rgx.match(position_seg)
            if position_match:
                position_match_dct = position_match.groupdict()
                position = int(position_match_dct['position'])
                # position_type = position_match_dct.get('position_type', '')
                position_info['positions'].append(position)
                position_info['positions_type'].append(position_seg.strip(','))

        return position_info

    def get_elem_info(self):
        sum_mod_elem_dct = {}
        for mod in self.mod_info:
            mod_cv = mod["cv"]
            mod_count = mod['count']
            mod_elem_dct = cv_elements_info[mod_cv]
            for elem in mod_elem_dct:
                if elem in sum_mod_elem_dct:
                    sum_mod_elem_dct[elem] += mod_count * mod_elem_dct[elem]
                else:
                    sum_mod_elem_dct[elem] = mod_count * mod_elem_dct[elem]

        return sum_mod_elem_dct

    def to_json(self):
        mod_json_str = json.dumps(self.mod_info)
        print(mod_json_str)
        is_valid = self.validator.is_valid(json.loads(mod_json_str))
        if is_valid:
            logger.debug(f"Schema check PASSED.")
            return mod_json_str
        else:
            raise ValueError(f"Schema check FAILED.")

    def to_mass_shift(self, angle_brackets: bool = True) -> str:
        if self.mod_level >= 0:
            sum_mod_elem_dct = self.get_elem_info()
            delta = 0
            for elem in sum_mod_elem_dct:
                if elem in elem_nominal_info:
                    delta += elem_nominal_info[elem] * sum_mod_elem_dct[elem]

            mod_str = f'{delta:+}'
            if angle_brackets:
                mod_str = self.__add_angle_brackets__(mod_str)
            return mod_str
        else:
            raise ValueError(f'Modification level not fit. required level >= 0, received: {self.mod_level}')

    def to_elem_shift(self, angle_brackets: bool = True) -> str:
        if self.mod_level >= 1:
            sum_mod_elem_dct = self.get_elem_info()
            mod_str_lst = []
            mod_elem_lst = ["O", "N", "S", "H", 'Na']
            for mod_elem in mod_elem_lst:
                if mod_elem in sum_mod_elem_dct:
                    mod_str_lst.append(f'{sum_mod_elem_dct[mod_elem]:+}{mod_elem}')

            mod_str = f'{",".join(mod_str_lst)}'
            if angle_brackets:
                mod_str = self.__add_angle_brackets__(mod_str)
            return mod_str
        else:
            raise ValueError(f'Modification level not fit. required level >= 1, received: {self.mod_level}')

    def to_mod_count(self, angle_brackets: bool = True) -> str:
        if self.mod_level >= 2:
            mod_output_lst = []
            for mod in self.mod_info:
                if mod["cv"] != 'DB':
                    mod_output_lst.append(f'{mod["count"]}{mod["cv"]}')
            mod_str = f'{",".join(mod_output_lst)}'
            if angle_brackets:
                mod_str = self.__add_angle_brackets__(mod_str)
            return mod_str
        else:
            raise ValueError(f'Modification level not fit. required level >= 2, received: {self.mod_level}')

    def to_mod_position(self, angle_brackets: bool = True) -> str:
        if self.mod_level >= 3:
            mod_output_lst = []
            for mod in self.mod_info:
                positions_lst = [str(i) for i in mod["positions"]]
                positions = '{' + ','.join(positions_lst) + '}'
                if mod["cv"] == 'DB':
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
            raise ValueError(f'Modification level not fit. required level >= 3, received: {self.mod_level}')

    def to_mod_position_type(self, angle_brackets: bool = True) -> str:
        if self.mod_level >= 4:
            mod_output_lst = []
            for mod in self.mod_info:
                positions = '{' + ','.join(mod["positions_type"]) + '}'
                if mod["cv"] == 'DB':
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
            raise ValueError(f'Modification level not fit. required level >= 4, received: {self.mod_level}')

    def to_all_levels(self, angle_brackets: bool = True) -> list:

        mod_output_lst = []
        if self.mod_level >= 0:
            mod_output_lst.append(mod_obj.to_mass_shift(angle_brackets=angle_brackets))
        if self.mod_level >= 1:
            mod_output_lst.append(mod_obj.to_elem_shift(angle_brackets=angle_brackets))
        if self.mod_level >= 2:
            mod_output_lst.append(mod_obj.to_mod_count(angle_brackets=angle_brackets))
        if self.mod_level >= 3:
            mod_output_lst.append(mod_obj.to_mod_position(angle_brackets=angle_brackets))
        if self.mod_level >= 4:
            mod_output_lst.append(mod_obj.to_mod_position_type(angle_brackets=angle_brackets))
        return mod_output_lst


if __name__ == '__main__':

    mod_code_lst = [
        r'<+46>',
        r'<+3O,-2H>',
        r'<2OH,Ke>',
        r'<2OH{8R,11S},Ke{14}>',
        r'<{5,9,12,15},2OH{8,11},Ke{14}>',
        r'<{5,9,12,15},2OH{8R,11S},Ke{14}>',
        r'<{5Z,9E,12E,15E},2OH{8,11},Ke{14}>',
        r'<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>',
    ]

    for usr_mod_code in mod_code_lst:
        print(usr_mod_code)

        mod_obj = Modifications(usr_mod_code)

        print(mod_obj.to_json())

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
        print('to all levels', mod_obj.to_all_levels())
        print('to all levels', mod_obj.to_all_levels(angle_brackets=False))
