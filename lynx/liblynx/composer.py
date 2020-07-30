# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
#
# LipidLynxX is using GPL V3 License
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: a data transfer hub to support integration of large scale lipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from lynx.controllers.converter import Converter
from lynx.controllers.decoder import Decoder
from lynx.controllers.encoder import Encoder
from lynx.models.lipid import (
    LipidClassType,
    ModType,
    ResidueType,
    LipidType,
    SiteType,
    DBType,
)


def construct_lipid(name):
    decoder = Decoder()
    extracted_info = decoder.extract(name)
    if extracted_info:
        for p in extracted_info:
            p_info = extracted_info[p]
            for in_r in p_info:
                r_info = p_info[in_r]
                segments = r_info.get("SEGMENTS", {})
                l_main_class = segments.get("CLASS", [""])[0]
                lipid_class = LipidClassType(main_class=l_main_class)
                residues_lst = []

                print(lipid_class)
                print(lipid_class.main_class)
                print(lipid_class.lmsd_sub_class)
                print(lipid_class.lmsd_sub_class)


if __name__ == "__main__":
    usr_lipid = "PC(16:0/20:4(7OH,9OH))"
    construct_lipid(usr_lipid)
