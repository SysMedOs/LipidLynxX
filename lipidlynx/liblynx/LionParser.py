# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from io import BytesIO
import base64

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem, rdMolDescriptors

from lipidlynx.controllers.Logger import logger
from lipidlynx.liblynx.LipidNomenclature import ParserFA, ParserPL
from lipidlynx.models.DefaultParams import default_cfg_path
from lipidlynx.liblynx.Converter import Converter


def parse_epilion(abbr: str) -> dict:

    fa_decoder = ParserFA()
    pl_decoder = ParserPL()

    info_dct = {}

    converter = Converter(default_cfg_path)
    epilion_id = converter.convert_abbr(abbr)

    if fa_decoder.is_fa(epilion_id):
        smi = fa_decoder.get_smi_fa(epilion_id)
        logger.info(epilion_id + ": " + smi)
    elif pl_decoder.is_pl(epilion_id):
        smi = pl_decoder.get_smi_pl(epilion_id)
        logger.info(epilion_id + ": " + smi)
    else:
        logger.info(f"Can NOT parse abbreviation: {epilion_id}")

    try:
        mol = Chem.MolFromSmiles(smi)
        AllChem.Compute2DCoords(mol)
        # m_mass = Descriptors.MolWt(mol)
        m_exact_mass = rdMolDescriptors.CalcExactMolWt(mol)
        m_formula = rdMolDescriptors.CalcMolFormula(mol)
        img = Draw.MolToImage(mol, size=(600, 400))
        img_io = BytesIO()
        img.save(img_io, format="png")
        img_io.seek(0)
        img.save(img_io, format="png")
        img_data = base64.b64encode(img_io.getbuffer())
        img_data_url = r"data:image/png;base64," + img_data.decode("utf-8")

        info_dct["id"] = epilion_id
        info_dct["formula"] = m_formula
        info_dct["exact_mass"] = "%.4f" % m_exact_mass
        info_dct["img"] = img_data_url

    except Exception as e:
        logger.error(f"! FAILED: {epilion_id}")
        logger.error(f"! FAILED to generate structure from SMILES: {smi}")
        logger.error(e)

    return info_dct
