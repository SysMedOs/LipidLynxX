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

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from epilion.libLION.DefaultParams import logger, abbr_cfg_df
from epilion.libLION.LipidNomenclature import ParserFA, ParserPL
from epilion.libLION.AbbrElemCalc import ElemCalc
from epilion.libLION.Converter import Converter

elem_calc = ElemCalc()
converter = Converter(abbr_df=abbr_cfg_df)


@dataclass
class Lipid(object):

    _id: str = field(init=False, repr=False)
    _lipidclass: str = field(init=False, repr=False)
    _formula: str = field(init=False, repr=False)
    _elements: dict = field(init=False, repr=False)
    _exactmass: float = field(init=False, repr=False)
    _smiles: str = field(init=False, repr=False)
    _properties: dict = field(init=False, repr=False)
    _spectra: dict = field(init=False, repr=False)

    input_abbr: str = None

    def __post_init__(self):
        self.__setattr__('_id', converter.convert_abbr(self.input_abbr))
        self.__setattr__('_lipidclass', self.id.split('(')[0])
        self.__setattr__('_formula', elem_calc.get_formula(self.id)[0])
        self.__setattr__('_elements', elem_calc.get_formula(self.id)[1])
        self.__setattr__('_exactmass', elem_calc.get_exactmass(self.elements))
        self.smiles = ''
        self._properties = {}
        self._spectra = {}

        self.__slots__ = ['id', 'lipidclass', 'formula', 'elements', 'exactmass', 'smiles', 'properties', 'spectra']

        self.input_abbr = self.input_abbr

    @property
    def id(self) -> str:
        return self._id

    @property
    def lipidclass(self):
        return self._lipidclass

    @property
    def elements(self):
        return self._elements

    @property
    def formula(self):
        return self._formula

    @property
    def exactmass(self):
        return self._exactmass

    @property
    def smiles(self) -> str:
        return self._smiles

    @smiles.setter
    def smiles(self, smi_str: str):
        self._smiles = smi_str

    @property
    def properties(self) -> dict:
        return self._properties

    @properties.setter
    def properties(self, properties_dct: dict):
        self._properties = properties_dct

    @property
    def spectra(self) -> dict:
        return self._spectra

    @spectra.setter
    def spectra(self, spectra_dct: dict):
        self._spectra = spectra_dct


@dataclass
class FA(Lipid):

    input_abbr: str = None

    def __post_init__(self):
        super().__post_init__()


@dataclass
class LCB(Lipid):

    input_abbr: str = None

    def __post_init__(self):
        super().__post_init__()


@dataclass
class PL(Lipid):
    _headgroup: str = field(init=False, repr=False)
    _bulk: str = field(init=False, repr=False)
    _fa1: FA = field(init=False, repr=False)
    _fa2: FA = field(init=False, repr=False)

    input_abbr: str = None

    def __post_init__(self):
        super().__post_init__()
        _hg = self.lipidclass.strip('Lyso')
        self.__setattr__('_headgroup', _hg.strip('L'))
        self.__setattr__('_bulk', 'PC(34:2)')
        self.__setattr__('_fa1', FA('FA16:0'))
        self.__setattr__('_fa2', FA('FA18:2'))

    @property
    def headgroup(self):
        return self._headgroup

    @property
    def bulk(self):
        return self._bulk

    @property
    def fa1(self) -> FA:
        return self._fa1

    @property
    def fa2(self) -> FA:
        return self._fa2


@dataclass
class SP(Lipid):
    _headgroup: str = field(init=False, repr=False)
    _bulk: str = field(init=False, repr=False)
    _lcb: LCB = field(init=False, repr=False)
    _fa: FA = field(init=False, repr=False)

    input_abbr: str = None

    def __post_init__(self):
        super().__post_init__()
        self.__setattr__('_headgroup', self.lipidclass)
        self.__setattr__('_bulk', 'SP(d-34:2)')
        self.__setattr__('_lcb', LCB('d-16:0'))
        self.__setattr__('_fa', FA('FA18:2'))

    @property
    def headgroup(self):
        return self._headgroup

    @property
    def bulk(self):
        return self._bulk

    @property
    def lcb(self) -> LCB:
        return self._lcb

    @property
    def fa(self) -> FA:
        return self._fa


@dataclass
class GL(Lipid):

    _bulk: str = field(init=False, repr=False)
    _fa1: FA = field(init=False, repr=False)
    _fa2: FA = field(init=False, repr=False)
    _fa3: FA = field(init=False, repr=False)

    input_abbr: str = ''

    def __post_init__(self):
        super().__post_init__()
        self.__setattr__('_bulk', 'PC(34:2)')
        self.__setattr__('_fa1', FA('FA16:0'))
        self.__setattr__('_fa2', FA('FA18:2'))
        self.__setattr__('_fa3', FA('FA18:2'))

    @property
    def bulk(self):
        return self._bulk

    @property
    def fa1(self) -> FA:
        return self._fa1

    @property
    def fa2(self) -> FA:
        return self._fa2

    @property
    def fa3(self) -> FA:
        return self._fa3


class LipidEncoder(json.JSONEncoder):
    def default(self, obj):
        if is_dataclass(obj):
            return asdict(obj)
        return super().default(obj)


class LipidEncoderLite(json.JSONEncoder):
    def default(self, obj):
        if is_dataclass(obj):
            in_dct = asdict(obj)
            out_dct = copy.deepcopy(in_dct)
            for l in in_dct:
                if is_dataclass(out_dct[l]):
                    in_dct[l] = asdict(in_dct[l])
                    out_dct[l] = copy.deepcopy(in_dct[l])

                if isinstance(out_dct[l], dict):
                    if in_dct[l]:
                        for m in in_dct[l]:
                            if m.strip('_') in ['smiles', 'properties', 'spectra']:
                                if not in_dct[l][m] or in_dct[l][m] == '':
                                    del out_dct[l][m]
                                if isinstance(in_dct[l][m], dict):
                                    if not not in_dct[l][m]:
                                        del out_dct[l][m]
                            elif m.strip('_') in ['input_abbr']:
                                del out_dct[l][m]
                    else:
                        print(f'del {l}')
                        del out_dct[l]
                else:
                    if l.strip('_') in ['smiles', 'properties', 'spectra']:
                        if not in_dct[l] or in_dct[l] == '':
                            del out_dct[l]
                    elif l.strip('_') in ['input_abbr']:
                        del out_dct[l]

            return out_dct
        return super().default(obj)


if __name__ == '__main__':

    usr_pl = r'PC(16:0_18:2)'
    usr_tg = r'TG(14:0_16:0_18:2)'
    pl = PL(input_abbr=usr_pl)
    print(pl.input_abbr)
    pl.input_abbr = 'x'
    pl.smiles = r'CCCCCCCCCCCCCCCCCCCCCCO'
    print(pl.input_abbr)
    print(pl.formula)
    print(pl.elements)
    print(pl.smiles)
    print(pl.fa1)
    print(type(pl.fa1), isinstance(pl.fa1, FA), isinstance(pl.fa1, Lipid))
    fa1 = pl.fa1
    print(type(pl.fa1.elements), pl.fa1.elements)
    pl.charge = 'neg'
    pl.bug = 'BUG'
    logger.info(asdict(pl))
    print(pl.bug)
    tg = GL(input_abbr=usr_tg)
    logger.info(asdict(tg))

    print(asdict(tg).get('smi', 'No SMILES'), type(tg.smiles))

    logger.info(json.dumps(pl, cls=LipidEncoder))
    logger.info(json.dumps(pl, cls=LipidEncoderLite))

    logger.info(json.dumps(tg, cls=LipidEncoder))
    logger.info(json.dumps(tg, cls=LipidEncoderLite))

