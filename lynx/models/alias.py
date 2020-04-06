# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json

import regex as re

from lynx.models.defaults import logger, default_alias_file
from lynx.utils.file_readers import get_json


class Alias(object):

    def __init__(self, alias_file: str = default_alias_file):
        self.alias_file = alias_file
        self._raw_cv = {}
        self._raw_cv = self.__load_raw__()
        self._info = self.__load__()

    def __load__(self) -> dict:
        alias_info = {}
        raw_alias_info = self.__load_raw__()
        for alias_category in raw_alias_info:
            raw_alias_cat_info = raw_alias_info[alias_category]
            alias_cat_info = {}
            for alias_class in raw_alias_cat_info:
                raw_alias_class_info = raw_alias_cat_info[alias_class]
                alias_class_info = {}
                for abbr in raw_alias_class_info:
                    alias_lst = raw_alias_class_info[abbr]
                    for alias in alias_lst:
                        alias_class_info[re.compile(alias)] = abbr
                alias_cat_info[alias_class] = alias_class_info
            alias_info[alias_category] = alias_cat_info
        return alias_info

    def __load_raw__(self) -> dict:
        if not self._raw_cv:
            return get_json(self.alias_file)
        else:
            return self._raw_cv

    # def __str__(self):
    #     return json.dumps(self.info)
    #
    # def __repr__(self):
    #     return self.__str__

    @property
    def info(self) -> dict:
        return self._info

    @property
    def raw_cv(self) -> dict:
        return self._raw_cv


if __name__ == '__main__':

    usr_alias = r'lynx/configurations/defined_alias.json'

    alias_obj = Alias(usr_alias)
    print(alias_obj.info)
    print(alias_obj)
