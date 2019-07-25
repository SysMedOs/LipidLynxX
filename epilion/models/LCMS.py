# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import copy
import json
import os.path
from dataclasses import dataclass, asdict, field, is_dataclass
from typing import Dict, Tuple, List

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from epilion.libLION.DefaultParams import logger, abbr_cfg_df
from epilion.libLION.LipidNomenclature import ParserFA, ParserPL
from epilion.libLION.AbbrElemCalc import ElemCalc
from epilion.libLION.Converter import Converter


@dataclass
class Spectrum:
    mz: List[float] = field(default_factory=list)
    i: List[float] = field(default_factory=list)
