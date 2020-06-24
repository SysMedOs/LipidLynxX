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
import re
from typing import Dict, List

from pydantic import BaseModel, constr


lipid_name_rgx_str = r"^\s*.{2,512}\s*$"
LipidNameType = constr(regex=lipid_name_rgx_str)
level_rgx_str = r"^\s*(?P<level>[Bb][0-3]?|^[DSds]([0-5](.[0-3])?)?|^[Mm][Aa][Xx])\s*$"
level_rgx = re.compile(level_rgx_str)


class LvType(str):
    """
    Define LipidLynxX level as specific Type.
    """

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        if not isinstance(v, str):
            raise TypeError("string required")
        m = level_rgx.match(v.upper())
        if not m:
            raise ValueError("invalid LipidLynxX level format")
        else:
            lv = m.groupdict().get("level")
        return cls(lv)

    def __repr__(self):
        return f"LvType({super().__repr__()})"

    class Config:
        schema_extra = {"example": "B1"}


class LevelsType(BaseModel):
    """
    Define LipidLynxX level as specific Type.
    """

    levels: List[LvType]

    class Config:
        schema_extra = {"example": ["B1", "D1"]}


class LevelsData(BaseModel):
    """
    Define LipidLynxX levels for input.
    """

    levels: List[str]

    class Config:
        schema_extra = {"example": {"levels": ["B1", "D1"]}}


class FileType(str, Enum):
    xlsx = "xlsx"
    csv = "csv"


class StyleType(str, Enum):
    lipidlynxx = "LipidLynxX"
    comp_db = "COMP_DB"
    shorthand = "ShorthandNotation"


class InputStrData(BaseModel):
    data: LipidNameType

    class Config:
        schema_extra = {"example": {"data": "PLPC"}}


class InputListData(BaseModel):
    data: List[LipidNameType]

    class Config:
        schema_extra = {
            "example": {"data": ["DHA", "PLPC", "PI 16:0-20:4", "UNKNOWN_LIPID_1"]}
        }


class InputDictData(BaseModel):
    data: Dict[str, List[LipidNameType]]

    class Config:
        schema_extra = {
            "example": {
                "data": {
                    "Source_1": ["DHA", "PLPC", "PI 16:0-20:4", "UNKNOWN_LIPID_1"],
                    "Source_2": [
                        "FA20:4",
                        "PC 16:0_18:2",
                        "PI(16:0/20:4)",
                        "UNKNOWN_LIPID_2",
                    ],
                }
            }
        }


class ConvertedStrData(BaseModel):
    input: str
    output: str
    converted: List[str]
    skipped: str

    class Config:
        schema_extra = {
            "example": {
                "input": "PLPC",
                "output": "PC(16:0/18:2)",
                "converted": [["PLPC", "PC(16:0/18:2)"],],
                "skipped": "",
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
                "input": ["DHA", "PLPC", "PI 16:0-20:4", "UNKNOWN_LIPID"],
                "output": ["FA22:6", "PC(16:0/18:2)", "PI(16:0_20:4)"],
                "converted": [
                    ["DHA", "FA22:6"],
                    ["PLPC", "PC(16:0/18:2)"],
                    ["PI 16:0-20:4", "PI(16:0_20:4)"],
                ],
                "skipped": ["UNKNOWN_LIPID"],
            }
        }


class ConverterExportData(BaseModel):
    data: Dict[str, ConvertedListData]

    class Config:
        schema_extra = {
            "example": {
                "data": {
                    "Source_1": {
                        "input": ["DHA", "PLPC", "PI 16:0-20:4", "UNKNOWN_LIPID_1"],
                        "output": ["FA22:6", "PC(34:2)", "PI(36:4)"],
                        "converted": [
                            ["DHA", "FA22:6"],
                            ["PLPC", "PC(34:2)"],
                            ["PI 16:0-20:4", "PI(36:4)"],
                        ],
                        "skipped": ["UNKNOWN_LIPID_1"],
                    },
                    "Source_2": {
                        "input": [
                            "FA20:4",
                            "PC 16:0_18:2",
                            "PI(16:0/20:4)",
                            "UNKNOWN_LIPID_2",
                        ],
                        "output": ["FA20:4", "PC(34:2)", "PI(36:4)"],
                        "converted": [
                            ["FA20:4", "FA20:4"],
                            ["PC 16:0_18:2", "PC(34:2)"],
                            ["PI(16:0/20:4)", "PI(36:4)"],
                        ],
                        "skipped": ["UNKNOWN_LIPID_2"],
                    },
                }
            }
        }


class EqualizedLevelData(BaseModel):
    matched: Dict[str, dict]
    unmatched: Dict[str, dict]

    class Config:
        schema_extra = {"example": {"matched": {}, "unmatched": {},}}


class EqualizedData(BaseModel):
    equalized: Dict[str, EqualizedLevelData]
    skipped: Dict[str, List[str]]

    class Config:
        schema_extra = {
            "example": {
                "equalized": {"B1": {"matched": {}, "unmatched": {},},},
                "skipped": {"Source01": ["bad_ID", "Unknown_id"]},
            }
        }


class EqualizerExportData(BaseModel):
    data: EqualizedData

    class Config:
        schema_extra = {
            "example": {
                "data": {"equalized": {"B1": None,}, "skipped": ["UNKNOWN_LIPID_1"],},
            }
        }
