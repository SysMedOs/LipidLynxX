# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
#
# LipidLynxX is using GPL V3 License
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: a data transfer hub to support integration of large scale lipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json

import regex as re

from lynx.models.defaults import default_cv_file
from lynx.utils.file_handler import get_json
from lynx.utils.log import app_logger


class CV(object):
    def __init__(self, cv_file: str = default_cv_file, logger=app_logger):
        self.cv_file = cv_file
        self._info = self.__load__()
        self._raw_cv = self.__load_raw__()
        self.logger = logger

    def __load__(self) -> dict:
        cv_info = {}
        raw_cv_info = get_json(self.cv_file)
        for raw_cv in raw_cv_info:
            order = raw_cv.get("order", None)
            cv = raw_cv.get("cv", None)
            level = raw_cv.get("level", 0)
            elements_dct = raw_cv.get("elements", {})
            alias_lst = raw_cv.get("alias", [])
            if cv:
                if cv not in alias_lst:
                    alias_lst.append(cv)
                    alias_lst = list(set(alias_lst))
                for alias in alias_lst:
                    cv_info[alias] = {
                        "MATCH": re.compile("^" + alias + "$"),
                        "ORDER": order,
                        "CV": cv,
                        "LEVEL": level,
                        "ELEMENTS": elements_dct,
                    }

        return cv_info

    def __load_raw__(self) -> dict:
        raw_cv_info = {}
        raw_cv_dct = get_json(self.cv_file)
        for raw_cv in raw_cv_dct:
            order = raw_cv.get("order", None)
            cv = raw_cv.get("cv", None)
            level = raw_cv.get("level", 0)
            elements_dct = raw_cv.get("elements", {})
            alias_lst = raw_cv.get("alias", [])
            if cv:
                if cv not in alias_lst:
                    alias_lst.append(cv)
                    alias_lst = list(set(alias_lst))
                raw_cv_info[cv] = {
                    "ORDER": order,
                    "ALIAS": alias_lst,
                    "LEVEL": level,
                    "ELEMENTS": elements_dct,
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


if __name__ == "__main__":

    usr_cv = r"lynx/configurations/controlled_vocabularies.json"

    cv_obj = CV(usr_cv)
    print(cv_obj.info)
    print(cv_obj)
