import json
import pandas as pd
import re
from typing import Union, List

from jsonschema import Draft7Validator

from lipidlynx.models.defaults import lynx_schema, cv_elements_info, logger, elem_nominal_info
from lipidlynx.models.modification import Modifications
from lipidlynx.models.patterns import fa_rgx
from lipidlynx.controllers.GeneralFunctions import get_abs_path


class FattyAcid(object):

    def __init__(self, lipid_code: str, db: int = 0):

        self.lipid_code = lipid_code.strip('FA')

        self.fa_info_dct = self.__post_init__()

        self.mod_info = self.fa_info_dct.get('mod_obj', None)
        if db and db < self.fa_info_dct['db']:
            self.db_count = db
        else:
            self.db_count = self.fa_info_dct['db']

    def __post_init__(self):
        
        fa_info_dct = {}
        fa_match = fa_rgx.match(self.lipid_code)

        if fa_match:
            fa_info_dct = fa_match.groupdict()
            mod_code = fa_info_dct.get('mod', None)
            if mod_code:
                mod_obj = Modifications(mod_code)
                fa_info_dct.update({'mod_obj': mod_obj})
            return fa_info_dct
        else:
            raise ValueError(f'Can not parse FA sting: {self.lipid_code}')

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
        is_valid = self.validator.is_valid(json.loads(mod_json_str))
        if is_valid:
            logger.debug(f"Schema check PASSED.")
            return mod_json_str
        else:
            raise Exception(f"Schema check FAILED.")

    def to_mass_shift(self, angle_brackets: bool = True) -> str:
        sum_mod_elem_dct = self.get_elem_info()
        delta = 0
        for elem in sum_mod_elem_dct:
            if elem in elem_nominal_info:
                delta += elem_nominal_info[elem] * sum_mod_elem_dct[elem]

        mod_str = f'{delta:+}'
        if angle_brackets:
            mod_str = self.__add_angle_brackets__(mod_str)
        return mod_str

    def to_elem_shift(self, angle_brackets: bool = True) -> str:
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

    def to_mod_count(self, angle_brackets: bool = True) -> str:
        mod_output_lst = []
        for mod in self.mod_info:
            if mod["cv"] != 'DB':
                mod_output_lst.append(f'{mod["count"]}{mod["cv"]}')
        mod_str = f'{",".join(mod_output_lst)}'
        if angle_brackets:
            mod_str = self.__add_angle_brackets__(mod_str)
        return mod_str

    def to_mod_position(self, angle_brackets: bool = True) -> str:
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

    def to_mod_position_type(self, angle_brackets: bool = True) -> str:

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

    def to_all_levels(self, angle_brackets: bool = True) -> list:

        mod_output_lst = [
            self.to_mass_shift(angle_brackets=angle_brackets),
            self.to_elem_shift(angle_brackets=angle_brackets),
            self.to_mod_count(angle_brackets=angle_brackets),
            self.to_mod_position(angle_brackets=angle_brackets),
            self.to_mod_position_type(angle_brackets=angle_brackets)
        ]

        return mod_output_lst


if __name__ == '__main__':
    usr_mod_code = r'FA20:4<{5Z,9E,12E,15E},2OH{8R,11S},Ke{14}>'

    fa_obj = FattyAcid(usr_mod_code)

    # print(fa_obj.to_json())
    #
    # print('to mass shift', fa_obj.to_mass_shift())
    # print('to elem shift', fa_obj.to_elem_shift())
    # print('to mod_count', fa_obj.to_mod_count())
    # print('to mod_position', fa_obj.to_mod_position())
    # print('to mod_position_type', fa_obj.to_mod_position_type())
    # print('to mass shift', fa_obj.to_mass_shift(angle_brackets=False))
    # print('to elem shift', fa_obj.to_elem_shift(angle_brackets=False))
    # print('to mod_count', fa_obj.to_mod_count(angle_brackets=False))
    # print('to mod_position', fa_obj.to_mod_position(angle_brackets=False))
    # print('to mod_position_type', fa_obj.to_mod_position_type(angle_brackets=False))
    # print('to all levels', fa_obj.to_all_levels())
    # print('to all levels', fa_obj.to_all_levels(angle_brackets=False))
