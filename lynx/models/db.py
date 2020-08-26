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
    mod_db_level_lst,
    default_output_rules,
)
from lynx.utils.basics import get_abs_path
from lynx.utils.cfg_reader import api_version
from lynx.utils.log import app_logger
from lynx.utils.params_loader import load_output_rule
from lynx.utils.toolbox import check_json


class DB(object):
    def __init__(
        self,
        mod_info: dict,
        db: int = 0,
        schema: str = "lynx_mod",
        output_rules: dict = default_output_rules,
        nomenclature: str = "LipidLynxX",
        logger=app_logger,
    ):
        self.logger = logger
        self.nomenclature = nomenclature
        self.export_rule = load_output_rule(output_rules, nomenclature)
        self.mod_rule = self.export_rule.get("MODS", None)
        self.mod_rule_orders = self.mod_rule.get("MOD", {}).get("ORDER", [])
        self.mod_separators = self.export_rule.get("SEPARATORS", [])
        if not self.mod_rule:
            raise ValueError(
                f"Cannot find output rule for 'MODS' from nomenclature: {nomenclature}."
            )
        self.mod_info = mod_info.get("MOD_INFO", {})
        self.schema = schema
        self.type = "DB"
        self.db_level = str(mod_info.get("MOD_LEVEL", 0))
        with open(get_abs_path(lynx_schema_cfg[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(
                json.load(s_obj),
                resolver=RefResolver(
                    f"file://{core_schema_path}", referrer=core_schema
                ),
            )

        self.db_count = db
