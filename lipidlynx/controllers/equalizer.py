# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de


import re

import pandas as pd

from lipidlynx.controllers.general_functions import get_abs_path
from lipidlynx.models.log import logger
from lipidlynx.models.lipid import Lipid
from lipidlynx.models.residues import FattyAcid


class Equalizer(object):
    def __init__(self, file_path: str, level: str = "auto"):
        self.level = level
        abs_path = get_abs_path(file_path)
        if abs_path.lower().endswith(".xlsx"):
            self.df = pd.read_excel(abs_path)
        elif abs_path.lower().endswith(".csv"):
            self.df = pd.read_csv(abs_path)
        else:
            raise ValueError(f"Cannot read file {abs_path}")
        self.df.fillna("")
        self.header_lst = self.df.columns.to_list()

        self.data = self.df.to_dict(orient="list")

    def convert_col(self, col_name):
        id_lst = self.data.get(col_name, [])
        equalized_id_dct = {}
        skipped_id_lst = []
        if id_lst:
            for _id in id_lst:
                if isinstance(_id, str) and _id:
                    logger.info(f"Convert {_id} to level {self.level}")
                    logger.info(f"test {_id}")
                    _lipid = Lipid(lipid_code=_id)
                else:
                    _lipid = None
                if _lipid:
                    if _lipid.is_modified and self.level in _lipid.linked_ids:
                        equalized_id_dct[_lipid.linked_ids[self.level]] = _id
                    elif (
                        not _lipid.is_modified
                        and f"{self.level[0]}0" in _lipid.linked_ids
                    ):
                        equalized_id_dct[_lipid.linked_ids[f"{self.level[0]}0"]] = _id
                    else:
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
        for ref in equalized_dct:
            for _id in equalized_dct[ref]:
                if _id not in sum_equalized_dct:
                    sum_equalized_dct[_id] = {ref: equalized_dct[ref][_id]}
                else:
                    sum_equalized_dct[_id][ref] = equalized_dct[ref][_id]
        sum_matched_dct = {}
        no_match_dct = {}
        for _id in sum_equalized_dct:
            if len(sum_equalized_dct[_id]) > 1:
                sum_matched_dct[_id] = sum_equalized_dct[_id]
            else:
                no_match_dct[_id] = sum_equalized_dct[_id]

        output_dct = {
            "matched": sum_matched_dct,
            "equalized": no_match_dct,
            "skipped": converted_dct["skipped"],
        }

        return output_dct

    def export(self, out_path: str):
        output_dct = self.cross_match()
        xlsx_writer = pd.ExcelWriter(out_path)
        if "matched" in output_dct:
            matched_dct = output_dct["matched"]
            if matched_dct:
                out_matched_df = pd.DataFrame.from_dict(matched_dct, orient="index")
                out_matched_df.index.names = [f"ID@Lv_{self.level}"]
                out_matched_df.sort_index().to_excel(
                    xlsx_writer, sheet_name=f"matched_{self.level   }"
                )
        if "equalized" in output_dct:
            equalized_dct = output_dct["equalized"]
            if equalized_dct:
                pd.DataFrame.from_dict(
                    equalized_dct, orient="index"
                ).sort_index().to_excel(
                    xlsx_writer, sheet_name=f"unmatched_{self.level}"
                )
        if "skipped" in output_dct:
            skipped_dct = output_dct["skipped"]
            if skipped_dct:
                pd.DataFrame.from_dict(
                    skipped_dct, orient="index"
                ).T.sort_index().to_excel(
                    xlsx_writer, sheet_name="skipped", index=False
                )

        xlsx_writer.save()


if __name__ == "__main__":
    test_file = r"test/test_input/test_equalizer.csv"
    test_level_lst = [
        "B1",
        "B1.1",
        "B1.2",
        "B2",
        "B2.1",
        "B2.2",
        "D1",
        "D1.1",
        "D1.2",
        "D2",
        "D2.1",
        "D2.2",
        "D3",
        "D3.1",
        "D3.2",
    ]
    for lv in test_level_lst:
        equalizer = Equalizer(test_file, level=lv)
        equalizer.export(
            f"../../test/test_output/equalizer/test_level_mapper_output_{lv}.xlsx"
        )
        logger.info("FIN")
