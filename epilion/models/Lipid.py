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

from epilion.libLION.DefaultParams import logger, abbr_cfg_df
from epilion.libLION.LipidNomenclature import ParserFA, ParserPL
from epilion.libLION.AbbrElemCalc import ElemCalc
from epilion.libLION.Converter import Converter


@dataclass
class Lipid(object):

    _id: str = field(init=False, repr=False)
    _lipidclass: str = field(init=False, repr=False)
    _formula: str = field(init=False, repr=False)
    _formula_dct: dict = field(init=False, repr=False)
    _exactmass: float = field(init=False, repr=False)

    input_abbr: str = None
    elem_calc = ElemCalc()

    def __post_init__(self):
        self.id = self.input_abbr
        self.abbr = self.input_abbr
        self.formula = self.id
        self.formula_dct = self.id
        self.exactmass = self.id
        self.lipidclass = self.id
        __slots__ = ['id', 'lipidclass', 'formula', 'formula_dct', 'exactmass']

    @property
    def id(self) -> str:
        return self._id

    @id.setter
    def id(self, abbr: str):
        converter = Converter(abbr_df=abbr_cfg_df)
        self._id = converter.convert_abbr(abbr)

    @property
    def lipidclass(self):
        return self._lipidclass

    @lipidclass.setter
    def lipidclass(self, epilion_id):
        self._lipidclass = epilion_id.split('(')[0]

    @property
    def formula(self):
        return self._formula

    @formula.setter
    def formula(self, lion_abbr: str):
        self._formula = self.elem_calc.get_formula(lion_abbr)[0]

    @property
    def formula_dct(self):
        return self._formula_dct

    @formula_dct.setter
    def formula_dct(self, lion_abbr: str):
        self._formula_dct = self.elem_calc.get_formula(lion_abbr)[1]

    @property
    def exactmass(self):
        return self._exactmass

    @exactmass.setter
    def exactmass(self, lion_abbr: str):
        usr_elem_dct = self.elem_calc.get_neutral_elem(lion_abbr)
        self._exactmass = self.elem_calc.get_exactmass(usr_elem_dct)


@dataclass
class FA(Lipid):

    input_abbr: str = None

    def __post_init__(self):
        super().__post_init__()


@dataclass
class PL(Lipid):
    _headgroup: str = field(init=False, repr=False)

    input_abbr: str = None

    def __post_init__(self):
        super().__post_init__()
        self.headgroup = self._id

    @property
    def id(self) -> str:
        return self._id

    @id.setter
    def id(self, abbr: str):
        converter = Converter(abbr_df=abbr_cfg_df)
        self._id = converter.convert_abbr(abbr)

    @property
    def headgroup(self):
        return self._headgroup

    @headgroup.setter
    def headgroup(self, lipid_class: str):
        _hg = lipid_class.strip('Lyso')
        self._headgroup = _hg.strip('L')


@dataclass
class GL(Lipid):
    _headgroup: str = field(init=False, repr=False)

    input_abbr: str = ''

    def __post_init__(self):
        super().__post_init__()
        self.headgroup = self._id

    @property
    def headgroup(self):
        return self._headgroup

    @headgroup.setter
    def headgroup(self, lipid_class: str):
        _hg = lipid_class.strip('Lyso')
        self._headgroup = _hg.strip('L')


if __name__ == '__main__':

    usr_pl = r'PC(16:0_18:2)'
    usr_tg = r'TG(14:0_16:0_18:2)'
    pl = PL(input_abbr=usr_pl)
    print(pl.abbr)
    print(pl.input_abbr)
    pl.input_abbr = 'x'
    print(pl.input_abbr)
    print(pl.formula)
    pl.formula = 'y'
    print(pl.formula)
    pl.charge = 'neg'
    pl.bug = 'BUG'
    logger.info(asdict(pl))
    tg = GL(input_abbr=usr_tg)
    logger.info(asdict(tg))
