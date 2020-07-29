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

import regex as re

from lynx.models.defaults import default_alias_file
from lynx.utils.file_handler import get_json
from lynx.utils.log import app_logger


class Alias(object):
    def __init__(self, alias_file: str = default_alias_file, logger=app_logger):
        self.alias_file = alias_file
        self._raw_cv = {}
        self._raw_cv = self.__load_raw__()
        self._info = self.__load__()
        self._residue_alias = self._info.get("RESIDUE_ALIAS")
        self._lipid_alias = self._info.get("LIPID_ALIAS")
        self.logger = logger

    @staticmethod
    def __get_alias_info__(raw_alias_cat_info: dict) -> dict:
        alias_abbr_info = {}
        for alias_class in raw_alias_cat_info:
            raw_alias_class_info = raw_alias_cat_info[alias_class]
            for abbr in raw_alias_class_info:
                alias_lst = raw_alias_class_info[abbr]
                for alias in alias_lst:
                    # force add end fit to alias
                    if alias.endswith("$"):
                        pass
                    else:
                        alias += "\\s*$"
                    if re.match(r"[-_\dA-Z]{2,}", alias):  # if alias all uppercase
                        alias_abbr_info[re.compile(alias)] = abbr
                    else:
                        alias_abbr_info[re.compile(alias, re.IGNORECASE)] = abbr
        return alias_abbr_info

    def __load__(self) -> dict:
        alias_info = {}
        raw_alias_info = self.__load_raw__()
        residue_alias_info = {}
        lipid_alias_info = {}
        for alias_category in raw_alias_info:
            if alias_category.upper().startswith("RESIDUE"):
                raw_alias_cat_info = raw_alias_info[alias_category]
                residue_alias_info = self.__get_alias_info__(raw_alias_cat_info)
            elif alias_category.upper().startswith("LIPID"):
                raw_alias_cat_info = raw_alias_info[alias_category]
                lipid_alias_info = self.__get_alias_info__(raw_alias_cat_info)
            else:
                pass
        alias_info = {
            "RESIDUE_ALIAS": residue_alias_info,
            "LIPID_ALIAS": lipid_alias_info,
        }
        return alias_info

    def __load_raw__(self) -> dict:
        if not self._raw_cv:
            return get_json(self.alias_file)
        else:
            return self._raw_cv

    @property
    def info(self) -> dict:
        return self._info

    @property
    def residue_alias(self) -> dict:
        return self._residue_alias

    @property
    def lipid_alias(self) -> dict:
        return self._lipid_alias

    @property
    def raw_cv(self) -> dict:
        return self._raw_cv


if __name__ == "__main__":

    usr_alias = r"lynx/configurations/defined_alias.json"

    alias_obj = Alias(usr_alias)
    print(alias_obj.info)
    print(alias_obj)
