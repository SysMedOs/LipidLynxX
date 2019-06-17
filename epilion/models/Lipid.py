# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import os.path
from dataclasses import dataclass, asdict, field

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from epilion.libLION.DefaultParams import logger
from epilion.libLION.LipidNomenclature import ParserFA, ParserPL


@dataclass
class Lipid:

    _id: field(init=False, repr=False)
    _formula: str = field(init=False, repr=False)

    abbr: str = ''

    # def __init__(self, abbr):
    #     self.logger = logger
    #
    #     self.fa_decoder = ParserFA()
    #     self.pl_decoder = ParserPL()
    #
    #     # self.id = abbr
    #     self.formula = self.id
    #     self.exactmass = self.formula
    #     self.lipidclass = self.id

    def __post_init__(self):

        self.id = self.abbr
        self.formula = self.id
        self.exactmass = self.formula
        self.lipidclass = self.id

    @property
    def id(self) -> str:
        return self._id

    @id.setter
    def id(self, abbr: str):
        self._id = f'LION_{abbr}'

    @property
    def formula(self):
        return self._formula

    @formula.setter
    def formula(self, abbr):
        self._formula = f'Formula_{abbr}'

    # @property
    # def exactmass(self):
    #     return self.__exactmass
    #
    # @exactmass.setter
    # def exactmass(self, formula):
    #     self.__exactmass = f'Formula_{formula}'
    #
    # @property
    # def lipidclass(self):
    #     return self.__lipidclass
    #
    # @lipidclass.setter
    # def lipidclass(self, epilion_id):
    #     self.__lipidclass = epilion_id[:2]


if __name__ == '__main__':

    usr_abbr = 'Test'
    lipid = Lipid(usr_abbr)

    logger.info(lipid.abbr)
    logger.info(lipid.id)
    logger.info(lipid.formula)
    logger.info(asdict(lipid))
