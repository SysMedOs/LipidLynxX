# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from dataclasses import dataclass, field
from typing import List


@dataclass
class Spectrum:
    mz: List[float] = field(default_factory=list)
    i: List[float] = field(default_factory=list)
