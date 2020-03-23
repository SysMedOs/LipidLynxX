# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os

from jsonschema import Draft7Validator, RefResolver

from lynx.utils.toolbox import check_json

from .defaults import (
    api_version,
    mod_level_lst,
    hg_schema,
    hg_schema_path,
    fa_schema,
    fa_schema_path,
)


class Modifications(object):

    # __slots__ = ("info", "mods", "mods_site", "mods_info")

    def __init__(self, mods: dict):
        self.mods = mods

    def decode(self):
        for mod in self.mods:
            pass


class Residue(object):

    __slots__ = (
        "info",
        "link",
        "c_count",
        "db_count",
        "o_count",
        "is_modified",
        "sum_mods",
        "mods",
        "mods_site",
        "mods_info",
    )

    def __init__(self, residue_info: dict):
        self.info = residue_info

    @property
    def link(self) -> str:
        return self.info.get("LINK", "")

    @property
    def c_count(self) -> int:
        return int(self.info.get("NUM_C", "0"))

    @property
    def db_count(self) -> int:
        return int(self.info.get("DB", "0"))

    @property
    def o_count(self) -> int:
        return int(self.info.get("NUM_O", "0"))

    @property
    def is_modified(self) -> bool:
        is_mod = False
        sum_mods = self.info.get("SUM_MODS", [])
        mods = self.info.get("MOD", [])
        if sum_mods and mods:
            is_mod = True
        return is_mod

    @property
    def sum_mods(self) -> list:
        return self.info.get("SUM_MODS", [])

    @property
    def mods(self) -> object:
        return Modifications(self.info)

    @property
    def mods_site(self) -> list:
        return self.info.get("MOD_SITE", [])

    @property
    def mods_info(self) -> list:
        return self.info.get("MOD_INFO", [])


class LipidClass(object):
    def __init__(self):
        pass
