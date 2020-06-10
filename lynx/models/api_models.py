# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# LipidLynxX is Dual-licensed
#   For academic and non-commercial use: GPLv2 License:
#   For commercial use: please contact the SysMedOs team by email.
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: lipid annotations converter for large scale lipidomics and epilipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from enum import Enum
from typing import List

from pydantic import BaseModel


class FileType(str, Enum):
    xlsx = "xlsx"
    csv = "csv"


class StyleName(str, Enum):
    lipidlynxx = "LipidLynxX"
    comp_db = "COMP_DB"


class InputListData(BaseModel):
    lipid_names: List[str]

    class Config:
        schema_extra = {
            "example": {
                "lipid_names": ["DHA", "PLPC", "PI 16:0-20:4", "UNKOWN"]
            }
        }


class ExportListData(BaseModel):
    input: List[str]
    output: List[str]
    converted: List[List[str]]
    skipped: List[str]

    class Config:
        schema_extra = {
            "example": {
                "input": [
                    "DHA",
                    "PLPC",
                    "PI 16:0-20:4",
                    "UNKOWN"
                ],
                "output": [
                    "FA22:6",
                    "PC(16:0/18:2)",
                    "PI(16:0_20:4)"
                ],
                "converted": [
                    [
                        "DHA",
                        "FA22:6"
                    ],
                    [
                        "PLPC",
                        "PC(16:0/18:2)"
                    ],
                    [
                        "PI 16:0-20:4",
                        "PI(16:0_20:4)"
                    ]
                ],
                "skipped": ["UNKOWN"]
            }
        }
