# -*- coding: utf-8 -*-
#
# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json

import pandas as pd

from lipidlynx.models.log import logger
from lipidlynx.controllers.InitParams import (
    load_cfg_info,
    build_parser,
    build_mod_parser,
)
from lipidlynx.controllers.general_functions import get_abs_path

# Define default values across LipidLynx

api_version = "0.1"

# load default values from files defined in config.ini
# following parameters generated will be used as global values
with open(get_abs_path(r"lipidlynx/configurations/CV.json"), "r") as cv_js:
    cv_alias_json = json.load(cv_js)

cv_order_list = []
cv_alias_info = {}
cv_elements_info = {}

for _mod in cv_alias_json:
    cv_alias_info[_mod["cv"]] = _mod["alias"]
    cv_order_list.append(_mod["cv"])
    cv_elements_info[_mod["cv"]] = _mod["elements"]

default_cfg_path = "config.ini"
cfg_info_dct = load_cfg_info(cfg_path=default_cfg_path)
class_rgx_dct, rgx_class_dct = build_parser(cfg_info_dct["rules"])
cv_rgx_dct = build_mod_parser(cv_alias_info)
mod_cfg_df = pd.read_csv(cfg_info_dct["mod_cfg"], index_col=0, na_values=None)
abbr_cfg_df = pd.read_excel(cfg_info_dct["abbr_cfg"])

lipid_level_lst = ["B", "D", "S"]
mod_level_lst = ["0", "1", "2", "3", "4", "5"]
db_level_lst = [".0", "0.1", "0.2"]
mod_db_level_lst = [
    "0",
    "0.1",
    "0.2",
    "1",
    "1.1",
    "1.2",
    "2",
    "2.1",
    "2.2",
    "3",
    "3.1",
    "3.2",
    "4",
    "4.1",
    "4.2",
    "5",
    "5.1",
    "5.2",
]

lynx_schema_cfg = {
    "lynx_mod": r"lipidlynx/models/schema/lynx_mod.schema.json",
    "lynx_fa": r"lipidlynx/models/schema/lynx_fa.schema.json",
    "lynx_hg": r"lipidlynx/models/schema/lynx_hg.schema.json",
    "lynx_core": r"lipidlynx/models/schema/lynx_core.schema.json",
}

core_schema_path = get_abs_path(lynx_schema_cfg["lynx_core"])
with open(core_schema_path, "r") as core_obj:
    core_schema = json.load(core_obj)
hg_schema_path = get_abs_path(lynx_schema_cfg["lynx_hg"])
with open(hg_schema_path, "r") as hg_json_obj:
    hg_schema = json.load(hg_json_obj)
fa_schema_path = get_abs_path(lynx_schema_cfg["lynx_fa"])
with open(fa_schema_path, "r") as fa_json_obj:
    fa_schema = json.load(fa_json_obj)
mod_schema_path = get_abs_path(lynx_schema_cfg["lynx_mod"])
with open(mod_schema_path, "r") as mod_json_obj:
    mod_schema = json.load(mod_json_obj)

elem_nominal_info = {"H": 1, "N": 14, "O": 16, "S": 32, "Na": 23}

# # legacy params
# pa_hg_elem_dct = {"C": 0, "H": 3, "O": 4, "P": 1, "N": 0}
# pc_hg_elem_dct = {"C": 5, "H": 14, "O": 4, "P": 1, "N": 1}
# pe_hg_elem_dct = {"C": 2, "H": 8, "O": 4, "P": 1, "N": 1}
# pg_hg_elem_dct = {"C": 3, "H": 9, "O": 6, "P": 1, "N": 0}
# pi_hg_elem_dct = {"C": 6, "H": 13, "O": 9, "P": 1, "N": 0}
# pip_hg_elem_dct = {"C": 6, "H": 14, "O": 12, "P": 2, "N": 0}
# ps_hg_elem_dct = {"C": 3, "H": 8, "O": 6, "P": 1, "N": 1}
# tg_hg_elem_dct = {"C": 0, "H": 0, "O": 0, "P": 0, "N": 0}
# fa_hg_elem_dct = {"C": 0, "H": 0, "O": 0, "P": 0, "N": 0}
# glycerol_bone_elem_dct = {"C": 3, "H": 2}
# link_o_elem_dct = {"O": -1, "H": 2}
# link_p_elem_dct = {"O": -1}
#
# elem_shifts = {
#     "FA": fa_hg_elem_dct,
#     "O-": link_o_elem_dct,
#     "P-": link_p_elem_dct,
#     "PA": pa_hg_elem_dct,
#     "PC": pc_hg_elem_dct,
#     "PE": pe_hg_elem_dct,
#     "PG": pg_hg_elem_dct,
#     "PI": pi_hg_elem_dct,
#     "PS": ps_hg_elem_dct,
#     "PIP": pip_hg_elem_dct,
#     "TG": tg_hg_elem_dct,
#     "BB": glycerol_bone_elem_dct,  # BB for glycerol BackBone
# }
#
# elem_info = {
#     "H": [(1.0078250321, 0.999885), (2.0141017780, 0.0001157)],
#     "D": [(2.0141017780, 0.0001157)],
#     "C": [(12.0, 0.9893), (13.0033548378, 0.0107)],
#     "N": [(14.0030740052, 0.99632), (15.0001088984, 0.00368)],
#     "O": [(15.9949146221, 0.99757), (16.99913150, 0.00038), (17.9991604, 0.00205)],
#     "Na": [(22.98976967, 1.0)],
#     "P": [(30.97376151, 1.0)],
#     "S": [
#         (31.97207069, 0.9493),
#         (32.97145850, 0.0076),
#         (33.96786683, 0.0429),
#         (35.96708088, 0.0002),
#     ],
#     "K": [(38.9637069, 0.932581), (39.96399867, 0.000117), (40.96182597, 0.067302)],
# }
#
#
# pl_smi_info = {
#     "PA": r"OP(O)(OCC(",
#     "PC": r"[O-]P(OCC[N+](C)(C)C)(OCC(",
#     "PC-CH3": r"[O]P(OCC[N](C)C)(OCC(",
#     "PE": r"OP(OCCN)(OCC(",
#     "PG": r"OP(OCC(O)CO)(OCC(",
#     "PS": r"OP(OCC(C(O)=O)N)(OCC(",
#     "PI": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O))(OCC(",
#     "PIP": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(OCC(",
#     "PI4P": r"OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(OCC(",
#     "gly_part": r")C",
#     "pl_end": r")=O",
# }
#
# tg_smi_info = {"gly_part_a": r"[H]C(C", "gly_part_b": r")(", "gly_part_c": r")C"}
#
# lipid_class_alias_info = {
#     "O-a": {"CLASS": "O-", "RULE_CLASS": "FA"},
#     "O-p": {"CLASS": "P-", "RULE_CLASS": "FA"},
#     "cer": {"CLASS": "Cer", "RULE_CLASS": "Cer"},
#     "CER": {"CLASS": "Cer", "RULE_CLASS": "Cer"},
#     "GPA": {"CLASS": "PA", "RULE_CLASS": "PL"},
#     "GPCho": {"CLASS": "PC", "RULE_CLASS": "PL"},
#     "GPEtn": {"CLASS": "PE", "RULE_CLASS": "PL"},
#     "GPGro": {"CLASS": "PG", "RULE_CLASS": "PL"},
#     "GPIns": {"CLASS": "PI", "RULE_CLASS": "PL"},
#     "GPSer": {"CLASS": "PS", "RULE_CLASS": "PL"},
#     "PlsA": {"CLASS": "PA", "RULE_CLASS": "PL"},
#     "PlsCho": {"CLASS": "PC", "RULE_CLASS": "PL"},
#     "PlsEtn": {"CLASS": "PE", "RULE_CLASS": "PL"},
#     "PlsGro": {"CLASS": "PG", "RULE_CLASS": "PL"},
#     "PlsIns": {"CLASS": "PI", "RULE_CLASS": "PL"},
#     "PlsSer": {"CLASS": "PS", "RULE_CLASS": "PL"},
#     "GPC": {"CLASS": "PC", "RULE_CLASS": "PL"},
#     "GPE": {"CLASS": "PE", "RULE_CLASS": "PL"},
#     "GPG": {"CLASS": "PG", "RULE_CLASS": "PL"},
#     "GPI": {"CLASS": "PI", "RULE_CLASS": "PL"},
#     "GPS": {"CLASS": "PS", "RULE_CLASS": "PL"},
#     "MAG": {"CLASS": "MG", "RULE_CLASS": "GL"},
#     "DAG": {"CLASS": "DG", "RULE_CLASS": "GL"},
#     "TAG": {"CLASS": "TG", "RULE_CLASS": "GL"},
# }


logger.info("Default parameters loaded successfully.")
