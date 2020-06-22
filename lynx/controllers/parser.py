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

import re
from pydantic import parse_obj_as

from lynx.controllers.converter import Converter
from lynx.controllers.decoder import Decoder
from lynx.models.lipid import LipidType


def parse_lipid(lipid_name: str):
    lynx_converter = Converter()
    converted_results = lynx_converter.convert_str(input_str=lipid_name)
    converted_name = converted_results.output
    decoder = Decoder()
    extracted_info = decoder.extract(converted_name)
    parsed_info = {}
    if extracted_info:
        for p in extracted_info:
            p_info = extracted_info[p]
            for in_r in p_info:

                if re.search("#LipidLynxX", in_r, re.IGNORECASE):
                    parsed_info["id"] = converted_name
                    r_info = p_info[in_r]
                    segments = r_info.get("SEGMENTS", {})
                    parsed_info["lipid_class"] = segments.get("CLASS", [""])[0]
                    parsed_info["residues"] = r_info.get("RESIDUES", {})
                    mod_type_lst = segments.get("MOD_TYPE")
                    mod_type_lst = [mod for mod in mod_type_lst if mod]

                    exact_sn_position = segments.get("RESIDUE_SEPARATOR")
                    exact_sn_position_level = parsed_info["residues"].get(
                        "RESIDUES_SEPARATOR_LEVEL"
                    )
                    if exact_sn_position == ["/"] or exact_sn_position_level == "S":
                        parsed_info["exact_sn_position"] = True
                    else:
                        parsed_info["exact_sn_position"] = False
                    if mod_type_lst:
                        parsed_info["is_modified"] = True
                    else:
                        parsed_info["is_modified"] = False
    if parsed_info:
        export_info = parse_obj_as(LipidType, parsed_info)
    else:
        export_info = {"msg": f"failed to parse {lipid_name}"}

    return export_info
