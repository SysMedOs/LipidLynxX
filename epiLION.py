# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from LibLION.DefaultParams import logger
from LibLION.LipidNomenclature import ParserFA, ParserPL


def epilion2sdf(abbr_lst, save_sdf):

    fa_decoder = ParserFA()
    pl_decoder = ParserPL()

    info_dct = {}

    for abbr in abbr_lst:
        logger.info(abbr)
        if fa_decoder.is_fa(abbr):
            smi = fa_decoder.get_smi_fa(abbr)
            logger.info(abbr + ': ' + smi)
            info_dct[abbr] = smi
        elif pl_decoder.is_pl(abbr):
            smi = pl_decoder.get_smi_pl(abbr)
            logger.info(abbr + ': ' + smi)
            info_dct[abbr] = smi
        else:
            logger.info(f'Can NOT parse abbreviation: {abbr}')

    sdf_writer = Chem.SDWriter(open(save_sdf, mode='w'))

    for m in info_dct:
        smi = info_dct[m]
        try:
            mol = Chem.MolFromSmiles(smi)
            AllChem.Compute2DCoords(mol)
            mol.SetProp('_Name', m)
            _lpp_mass = Descriptors.MolWt(mol)
            _lpp_exactmass = rdMolDescriptors.CalcExactMolWt(mol)
            _lpp_formula = rdMolDescriptors.CalcMolFormula(mol)
            mol.SetProp('EXACT_MASS', '%.6f' % _lpp_exactmass)
            mol.SetProp('NOMINAL_MASS', '%.3f' % _lpp_mass)
            mol.SetProp('FORMULA', _lpp_formula)
            sdf_writer.write(mol)
        except Exception as e:
            logger.error(f'FAILED: {m}')
            logger.error(f'FAILED: {smi}')
            logger.error(e)


if __name__ == '__main__':

    fa_lst = [
        'FA18:0', '18:1', 'O-16:0', 'P-18:0',
        '20:4[4DB,2OH,1Ke]',
        '20:4[4DB{5,9,12,15},2OH{8,11},1Ke{14}]',
        '20:4[4DB{5Z,9E,12E,15E},2OH{8S,11R},1Ke{14}]',
        '20:4[4DB{5Z,9E,11Z,14Z},1OH{8S}]',
        # '9:0<CHO{@9C}>',
        # '20:1[PGA{8a,12b},1DB{13Z},1OH{15S}]'
    ]

    pl_lst = [
        r'PC(O-16:0/18:1)', r'PC(P-16:0_18:1)', r'PC(P-16:0/18:1)',
        'PC(16:0/20:4[4DB,2OH,1Ke])',
        'PC(16:0/20:4[4DB{5,9,12,15},2OH{8,11},1Ke{14}])',
        'PC(16:0/20:4[4DB{5Z,9E,12E,15E},2OH{8S,11R},1Ke{14}])',
    ]

    sum_lst = fa_lst + pl_lst

    epilion2sdf(sum_lst, 'temp/test.sdf')
    logger.info('FINISHED!')
