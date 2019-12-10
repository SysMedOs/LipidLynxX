# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from typing import Dict

from lipidlynx.models.log import logger


class LynxID:
    def __init__(self, lynx_id: str):
        self.lynx_id = lynx_id
        self.lynx_id_info = {}

    def add_ref_id(self, db: str, ref_id: str):
        ref_id_dct = self.lynx_id_info.get("resource_ids", {})
        if db in ref_id_dct:
            if ref_id_dct[db] != ref_id:
                logger.warning(
                    f"update ref_id in db: {db} from {ref_id_dct[db]} to {ref_id}"
                )
        ref_id_dct[db] = ref_id

    def add_ref_ids(self, db_ref_dct: Dict[str, str]):
        for db in db_ref_dct:
            self.add_ref_id(db, db_ref_dct[db])

    def create_json(self):
        pass
