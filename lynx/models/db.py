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
from typing import Dict, List, Union

from jsonschema import Draft7Validator, RefResolver
from natsort import natsorted
import regex as re

from lynx.models.defaults import (
    lynx_schema_cfg,
    core_schema,
    core_schema_path,
    default_output_rules,
    db_level_lst,
    mod_level_lst,
)
from lynx.utils.basics import get_abs_path
from lynx.utils.cfg_reader import api_version
from lynx.utils.log import app_logger
from lynx.utils.params_loader import load_output_rule
from lynx.utils.toolbox import check_json


class DB(object):
    def __init__(
        self,
        db_info: dict,
        schema: str = "lynx_db",
        output_rules: dict = default_output_rules,
        nomenclature: str = "LipidLynxX",
        logger=app_logger,
    ):
        self.logger = logger
        self.nomenclature = nomenclature
        self.export_rule = load_output_rule(output_rules, nomenclature)
        self.db_sites_rule = self.export_rule.get("DB_SITES", None)
        self.db_separators = self.export_rule.get("SEPARATORS", [])
        if not self.db_sites_rule:
            raise ValueError(
                f"Cannot find output rule for 'MODS' from nomenclature: {nomenclature}."
            )
        self.db_info = db_info.get("DB_INFO", {}).get("0.0_DB", {})
        self.schema = schema
        self.type = "DB"
        self.db_level = str(db_info.get("DB_LEVEL", 0))
        if self.db_level == "0":
            self.db_level = "0.0"
        with open(get_abs_path(lynx_schema_cfg[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(
                json.load(s_obj),
                resolver=RefResolver(
                    f"file://{core_schema_path}", referrer=core_schema
                ),
            )

        self.db_count = self.db_info.get("DB_COUNT", 0)
        self.db_site = self.to_db_site_list()
        self.db_site_info = self.to_db_site_info_list()
        self.sum_db_info = self.to_sum_info()

    def __str__(self):
        return self.to_json()

    def __repr__(self):
        return self.to_json()

    def to_db_site_list(self) -> List[str]:
        db_site_list = natsorted(self.db_info.get("DB_SITE", []))
        return db_site_list

    def to_db_site_info_list(self) -> List[str]:
        db_order = (
            self.export_rule.get("DB_SITES", {}).get("DB_INFO", {}).get("ORDER", {})
        )
        db_info_lst = []
        for db_seg in db_order:
            if db_seg in self.db_info:
                db_info_lst.append(self.db_info[db_seg])
        db_info_tp_lst = natsorted(list(zip(*db_info_lst)))
        db_info_seg_lst = ["".join(seg) for seg in db_info_tp_lst]
        return db_info_seg_lst

    def to_db_site(self) -> str:
        db_site_str = ",".join(self.db_site)
        return db_site_str

    def to_db_site_info(self) -> str:
        db_info_seg_lst = self.to_db_site_info_list()
        db_site_info_str = ",".join(db_info_seg_lst)
        return db_site_info_str

    def to_db_level(self, level: Union[int, float, str] = 0) -> str:
        db_str = ""
        if not isinstance(level, str):
            level = str(level)
        if float(level) > float(self.db_level):
            raise ValueError(
                f'Cannot convert to higher level than the db_level "{self.db_level}". Input:{level}'
            )

        if level in ["0", "0.0"]:
            db_str = ""
        elif level == "0.1":
            db_str = self.to_db_site()
        elif level == "0.2":
            db_str = self.to_db_site_info()
        else:
            raise ValueError(f"Currently not supported modification level: {level}")

        return db_str

    def to_all_levels(self, as_list: bool = False) -> Union[Dict[str, str], List[str]]:
        all_levels_dct = {}
        if self.db_level in db_level_lst:
            db_idx = db_level_lst.index(self.db_level)
            output_levels_lst = db_level_lst[: db_idx + 1]
        else:
            raise ValueError(f"DB level not supported: {self.db_level}")

        for level in output_levels_lst:
            all_levels_dct[level] = self.to_db_level(level)
        all_levels_info = all_levels_dct

        return all_levels_info

    def get_db_info(self) -> dict:
        db_js_dct = {
            "cv": "DB",
            "count": self.db_count,
            "site": [int(s) for s in self.db_site],
            "site_info": self.db_site_info,
            "order": 0,
            "elements": {},
            "mass_shift": 0,
        }
        return db_js_dct

    def to_sum_info(self):
        linked_ids = self.to_all_levels()
        db_id = linked_ids.get(self.db_level, "")
        if float(self.db_level) > 0 and db_id:
            sum_db_info_dct = {
                "api_version": api_version,
                "type": self.type,
                "id": db_id,
                "level": self.db_level,
                "linked_ids": linked_ids,
                "linked_levels": natsorted(list(linked_ids.keys())),
                "info": self.get_db_info(),
            }
        elif float(self.db_level) > 0 and not db_id:
            sum_db_info_dct = {
                "api_version": api_version,
                "type": self.type,
                "id": db_id,
                "level": self.db_level,
                "linked_ids": linked_ids,
                "linked_levels": natsorted(list(linked_ids.keys())),
                "info": self.get_db_info(),
            }
        elif float(self.db_level) == 0:
            sum_db_info_dct = {}
        else:
            raise ValueError(
                f"Cannot format DB code to level {self.db_level} "
                f"from input: {self.db_info}"
            )

        return sum_db_info_dct

    def to_json(self):
        db_json_str = json.dumps(self.sum_db_info)

        if check_json(
            validator=self.validator,
            json_obj=json.loads(db_json_str),
        ):
            return db_json_str
        else:
            raise Exception(f"Schema test FAILED. Schema {self.schema}")


if __name__ == "__main__":

    test_db_js = {
        "DB_LEVEL": 0.2,
        "DB_INFO": {
            "0.0_DB": {
                "DB_CV": "",
                "DB_LEVEL": 0.2,
                "DB_COUNT": 4,
                "DB_SITE": ["5", "8", "11", "14"],
                "DB_SITE_INFO": ["Z", "Z", "Z", "Z"],
                "DB_ORDER": 0.01,
            }
        },
    }

    db_obj = DB(db_info=test_db_js)
    print(db_obj.to_json())
