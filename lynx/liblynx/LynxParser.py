# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from io import BytesIO
import base64

# from rdkit import Chem
# from rdkit.Chem import Draw
# from rdkit.Chem import AllChem, rdMolDescriptors

from lynx.utils.log import logger
from ..liblynx.LipidNomenclature import ParserFA, ParserPL
from ..controllers.converter import Converter


def parse_lipidlynx(abbr: str) -> dict:

    fa_decoder = ParserFA()
    pl_decoder = ParserPL()

    info_dct = {}

    converted_dct = Converter.convert_string(abbr)
    lynx_id_lst = converted_dct.get("output", [])
    if len(lynx_id_lst) == 1:
        lynx_id = lynx_id_lst[0]
    else:
        lynx_id = ""

    if fa_decoder.is_fa(lynx_id):
        smi = fa_decoder.get_smi_fa(lynx_id)
        logger.info(lynx_id + ": " + smi)
    elif pl_decoder.is_pl(lynx_id):
        smi = pl_decoder.get_smi_pl(lynx_id)
        logger.info(lynx_id + ": " + smi)
    else:
        smi = ""
        logger.info(f"Can NOT parse abbreviation: {lynx_id}")

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

        info_dct["id"] = lynx_id
        info_dct["formula"] = m_formula
        info_dct["exact_mass"] = "%.4f" % m_exact_mass
        info_dct["img"] = img_data_url

    except Exception as e:
        logger.error(f"! FAILED: {lynx_id}")
        logger.error(f"! FAILED to generate structure from SMILES: {smi}")
        logger.error(e)

    return info_dct
