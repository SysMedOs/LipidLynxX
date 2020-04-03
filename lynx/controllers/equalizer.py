# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from typing import List, Union

import pandas as pd

from lynx.controllers.encoder import Encoder
from lynx.models.lipid import Lipid
from lynx.utils.file_readers import get_abs_path, create_equalizer_output
from lynx.utils.log import logger


class Equalizer(object):
    def __init__(self, input_data: Union[str, dict], level: Union[str, List[str]]):

        if isinstance(input_data, str):
            abs_path = get_abs_path(input_data)
            if abs_path.lower().endswith(".xlsx"):
                df = pd.read_excel(abs_path)
            elif abs_path.lower().endswith(".csv"):
                df = pd.read_csv(abs_path)
            else:
                raise ValueError(f"Cannot read file {abs_path}")
            df.fillna("")
            self.data = df.to_dict(orient="list")
        elif isinstance(input_data, dict):
            self.data = input_data
        else:
            raise ValueError(f"Not supported input {type(input_data)}")
        if isinstance(level, str):
            self.levels = [level]
        else:
            self.levels = level
        self.encoder = Encoder()
        self.header_lst = self.data.keys()

    def convert_col(self, col_name):
        id_lst = self.data.get(col_name, [])
        equalized_id_dct = {}
        skipped_id_lst = []
        if id_lst:
            for _id in id_lst:
                if isinstance(_id, str) and _id:
                    logger.info(f"Convert {_id} to level {self.levels}")
                    try:
                        _lynx_lv_id_dct = self.encoder.export_levels(_id, self.levels)
                        if _lynx_lv_id_dct:
                            equalized_id_dct[_id] = _lynx_lv_id_dct
                            logger.debug(
                                f"Converted input {_id} to level {self.levels}: {_lynx_lv_id_dct}"
                            )
                        else:
                            skipped_id_lst.append(_id)
                            logger.warning(
                                f"Cannot converted input {_id} to level {self.levels}"
                            )
                    except Exception as e:
                        logger.error(e)
                        skipped_id_lst.append(_id)
                else:
                    skipped_id_lst.append(_id)

        return {"equalized": equalized_id_dct, "skipped": skipped_id_lst}

    def convert_all(self):
        sum_dct = {"equalized": {}, "skipped": {}}
        for col_name in self.header_lst:
            _col_dct = self.convert_col(col_name)
            sum_dct["equalized"][col_name] = _col_dct["equalized"]
            sum_dct["skipped"][col_name] = _col_dct["skipped"]
        return sum_dct

    def cross_match(self):
        converted_dct = self.convert_all()
        equalized_dct = converted_dct["equalized"]
        sum_equalized_dct = {}
        no_match_dct = {}
        for lv in self.levels:
            lv_equalized_dct = {}
            for ref in equalized_dct:
                for source_id in equalized_dct[ref]:
                    lv_converted_id = equalized_dct[ref][source_id].get(lv, None)
                    if lv_converted_id:
                        if lv_converted_id not in lv_equalized_dct:
                            lv_equalized_dct[lv_converted_id] = {ref: source_id}
                        else:
                            lv_equalized_dct[lv_converted_id][ref] = source_id
            sum_matched_dct = {}
            for _id in lv_equalized_dct:
                if len(list(lv_equalized_dct[_id].keys())) > 1:
                    sum_matched_dct[_id] = lv_equalized_dct[_id]
                else:
                    no_match_dct[_id] = lv_equalized_dct[_id]
            if sum_matched_dct:
                sum_equalized_dct[lv] = sum_matched_dct
            else:
                pass

        output_dct = {}
        for m_lv in sum_equalized_dct:
            output_dct[f'Matched@Lv{m_lv}'] = sum_equalized_dct[m_lv]

        output_dct["Equalized"] = no_match_dct
        output_dct["Skipped"] = converted_dct["skipped"]

        return output_dct

    def export_dict(self):

        return self.cross_match()
