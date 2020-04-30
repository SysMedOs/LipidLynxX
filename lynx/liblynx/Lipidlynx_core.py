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

import os.path

# from rdkit import Chem
# from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from lynx.models.defaults import logger
from lynx.liblynx.LipidNomenclature import ParserFA, ParserPL


def lynx2sdf(abbr_lst, save_sdf):

    if isinstance(abbr_lst, str):
        try:
            if os.path.isfile(abbr_lst):
                logger.info(f"Try to open file: {abbr_lst}")
                with open(abbr_lst, "r") as infile_obj:
                    abbr_lst = infile_obj.readlines()
            else:
                logger.error(f"Can NOT __load__ input: {abbr_lst}")
                logger.info("!! END PROCESSING !!")
                exit()
        except Exception as e:
            logger.error(f"Can NOT __load__ input: {abbr_lst}")
            logger.error(e)

    fa_decoder = ParserFA()
    pl_decoder = ParserPL()

    info_dct = {}

    for abbr in abbr_lst:
        logger.info(abbr)
        if fa_decoder.is_fa(abbr):
            smi = fa_decoder.get_smi_fa(abbr)
            logger.info(abbr + ": " + smi)
            info_dct[abbr] = smi
        elif pl_decoder.is_pl(abbr):
            smi = pl_decoder.get_smi_pl(abbr)
            logger.info(abbr + ": " + smi)
            info_dct[abbr] = smi
        else:
            logger.info(f"Can NOT parse abbreviation: {abbr}")

    sdf_writer = Chem.SDWriter(open(save_sdf, mode="w"))

    for m in abbr_lst:
        if m in info_dct:
            smi = info_dct[m]
            try:
                mol = Chem.MolFromSmiles(smi)
                AllChem.Compute2DCoords(mol)
                mol.SetProp("_Name", m)
                m_mass = Descriptors.MolWt(mol)
                m_exactmass = rdMolDescriptors.CalcExactMolWt(mol)
                m_formula = rdMolDescriptors.CalcMolFormula(mol)
                mol.SetProp("EXACT_MASS", "%.6f" % m_exactmass)
                mol.SetProp("NOMINAL_MASS", "%.3f" % m_mass)
                mol.SetProp("FORMULA", m_formula)
                sdf_writer.write(mol)
            except Exception as e:
                logger.error(f"! FAILED: {m}")
                logger.error(f"! FAILED to generate structure from SMILES: {smi}")
                logger.error(e)
        else:
            logger.warning(f"!! Can NOT parse: {m}")


if __name__ == "__main__":

    test_file = r"../../test/test_input/test_names.txt"
    output_file = r"../../test/test_output/test_names_sdf.sdf"

    with open(test_file, "r") as input_obj:
        input_lst = input_obj.readlines()
        lynx2sdf(input_lst, output_file)

    lynx2sdf(test_file, output_file)
    lynx2sdf(input_lst, output_file)

    logger.info("FINISHED!")
