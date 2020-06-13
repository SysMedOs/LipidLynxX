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
from typing import Dict, List, Optional

from pydantic import BaseModel, Field, constr


lipid_name_rgx_str = r"^\s*.{2,512}\s*$"
LipidNameType = constr(regex=lipid_name_rgx_str)
level_rgx_str = r"^[Bb][0-3]?$|^[DSds]([0-5](.[0-3])?)?$|^[Mm][Aa][Xx]$"
LvType = constr(regex=level_rgx_str)


class FileType(str, Enum):
    xlsx = "xlsx"
    csv = "csv"


class StyleType(str, Enum):
    lipidlynxx = "LipidLynxX"
    comp_db = "COMP_DB"


class Levels(List[LvType]):
    class Config:
        schema_extra = {"example": ["B1", "D1"]}


class InputStrData(BaseModel):
    lipid_name: LipidNameType

    class Config:
        schema_extra = {"example": {"lipid_name": "PLPC"}}


class InputListData(BaseModel):
    lipid_names: List[LipidNameType]

    class Config:
        schema_extra = {
            "example": {"lipid_names": ["DHA", "PLPC", "PI 16:0-20:4", "Bad_Example_1"]}
        }


class InputDictData(BaseModel):
    data: Dict[str, List[LipidNameType]]

    class Config:
        schema_extra = {
            "example": {
                "data": {
                    "Source_1": ["DHA", "PLPC", "PI 16:0-20:4", "Bad_Example_1"],
                    "Source_2": [
                        "FA20:4",
                        "PC 16:0_18:2",
                        "PI(16:0/20:4)",
                        "Bad_Example_2",
                    ],
                }
            }
        }


class ConvertedListData(BaseModel):
    input: List[str]
    output: List[str]
    converted: List[List[str]]
    skipped: List[str]

    class Config:
        schema_extra = {
            "example": {
                "input": ["DHA", "PLPC", "PI 16:0-20:4", "UNKOWN"],
                "output": ["FA22:6", "PC(16:0/18:2)", "PI(16:0_20:4)"],
                "converted": [
                    ["DHA", "FA22:6"],
                    ["PLPC", "PC(16:0/18:2)"],
                    ["PI 16:0-20:4", "PI(16:0_20:4)"],
                ],
                "skipped": ["Bad_Example_1"],
            }
        }


class ConverterExportDictData(BaseModel):
    data: Dict[str, ConvertedListData]

    class Config:
        schema_extra = {
            "example": {
                "data": {
                    "Source_1": {
                        "input": ["DHA", "PLPC", "PI 16:0-20:4"],
                        "output": ["FA22:6", "PC(34:2)", "PI(36:4)"],
                        "converted": [
                            ["DHA", "FA22:6"],
                            ["PLPC", "PC(34:2)"],
                            ["PI 16:0-20:4", "PI(36:4)"],
                        ],
                        "skipped": ["Bad_Example_1"],
                    },
                    "Source_2": {
                        "input": ["FA20:4", "PC 16:0_18:2", "PI(16:0/20:4)"],
                        "output": ["FA20:4", "PC(34:2)", "PI(36:4)"],
                        "converted": [
                            ["FA20:4", "FA20:4"],
                            ["PC 16:0_18:2", "PC(34:2)"],
                            ["PI(16:0/20:4)", "PI(36:4)"],
                        ],
                        "skipped": ["Bad_Example_2"],
                    },
                }
            }
        }
