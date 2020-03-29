# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json

import regex as re

from lynx.models.defaults import default_cv_file
from lynx.utils.file_readers import get_json


class CV(object):

    def __init__(self, cv_file: str = default_cv_file):
        self.cv_file = cv_file
        self._info = self.__load__()
        self._raw_cv = self.__load_raw__()

    def __load__(self) -> dict:
        cv_info = {}
        raw_cv_info = get_json(self.cv_file)
        for raw_cv in raw_cv_info:
            order = raw_cv.get("order", None)
            cv = raw_cv.get('cv', None)
            level = raw_cv.get('level', 0)
            elements_dct = raw_cv.get('elements', {})
            alias_lst = raw_cv.get('alias', [])
            if cv:
                if cv not in alias_lst:
                    alias_lst.append(cv)
                    alias_lst = list(set(alias_lst))
                for alias in alias_lst:
                    cv_info[alias] = {
                        "MATCH": re.compile('^' + alias + '$'),
                        "ORDER": order,
                        "CV": cv,
                        "LEVEL": level,
                        "ELEMENTS": elements_dct
                    }

        return cv_info

    def __load_raw__(self) -> dict:
        raw_cv_info = {}
        raw_cv_dct = get_json(self.cv_file)
        for raw_cv in raw_cv_dct:
            order = raw_cv.get("order", None)
            cv = raw_cv.get('cv', None)
            level = raw_cv.get('level', 0)
            elements_dct = raw_cv.get('elements', {})
            alias_lst = raw_cv.get('alias', [])
            if cv:
                if cv not in alias_lst:
                    alias_lst.append(cv)
                    alias_lst = list(set(alias_lst))
                raw_cv_info[cv] = {
                    "ORDER": order,
                    "ALIAS": alias_lst,
                    "LEVEL": level,
                    "ELEMENTS": elements_dct
                }

        return raw_cv_info

    def __str__(self):
        return json.dumps(self.info)

    def __repr__(self):
        return self.__str__

    @property
    def info(self) -> dict:
        return self._info

    @property
    def raw_cv(self) -> dict:
        return self._raw_cv


if __name__ == '__main__':

    usr_cv = r'lynx/configurations/CV.json'

    cv_obj = CV(usr_cv)
    print(cv_obj.info)
    print(cv_obj)
