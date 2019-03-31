# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import logging

import pandas as pd

log_level = logging.DEBUG
logging.basicConfig(format='%(asctime)s-%(levelname)s - %(message)s', datefmt='%b-%d@%H:%M:%S', level=log_level)
logger = logging.getLogger('log')

# Define default values

mod_cfg_csv = r'../Configurations/Mod_cfg.csv'
mod_cfg_df = pd.read_csv(mod_cfg_csv, index_col=0, na_values=None)

logger.debug(mod_cfg_df)

pa_hg_elem_dct = {'C': 0, 'H': 3, 'O': 4, 'P': 1, 'N': 0}
pc_hg_elem_dct = {'C': 5, 'H': 14, 'O': 4, 'P': 1, 'N': 1}
pe_hg_elem_dct = {'C': 2, 'H': 8, 'O': 4, 'P': 1, 'N': 1}
pg_hg_elem_dct = {'C': 3, 'H': 9, 'O': 6, 'P': 1, 'N': 0}
pi_hg_elem_dct = {'C': 6, 'H': 13, 'O': 9, 'P': 1, 'N': 0}
pip_hg_elem_dct = {'C': 6, 'H': 14, 'O': 12, 'P': 2, 'N': 0}
ps_hg_elem_dct = {'C': 3, 'H': 8, 'O': 6, 'P': 1, 'N': 1}
tg_hg_elem_dct = {'C': 0, 'H': 0, 'O': 0, 'P': 0, 'N': 0}
fa_hg_elem_dct = {'C': 0, 'H': 0, 'O': 0, 'P': 0, 'N': 0}
glycerol_bone_elem_dct = {'C': 3, 'H': 2}
link_o_elem_dct = {'O': -1, 'H': 2}
link_p_elem_dct = {'O': -1}

elem_shifts = {'FA': fa_hg_elem_dct, 'O-': link_o_elem_dct, 'P-': link_p_elem_dct,
               'PA': pa_hg_elem_dct, 'PC': pc_hg_elem_dct, 'PE': pe_hg_elem_dct, 'PG': pg_hg_elem_dct,
               'PI': pi_hg_elem_dct, 'PS': ps_hg_elem_dct, 'PIP': pip_hg_elem_dct,
               'TG': tg_hg_elem_dct, 'BB': glycerol_bone_elem_dct,  # BB for glycerol BackBone
               }

elem_info = {'H': [(1.0078250321, 0.999885), (2.0141017780, 0.0001157)],
             'D': [(2.0141017780, 0.0001157)],
             'C': [(12.0, 0.9893), (13.0033548378, 0.0107)],
             'N': [(14.0030740052, 0.99632), (15.0001088984, 0.00368)],
             'O': [(15.9949146221, 0.99757), (16.99913150, 0.00038), (17.9991604, 0.00205)],
             'Na': [(22.98976967, 1.0)],
             'P': [(30.97376151, 1.0)],
             'S': [(31.97207069, 0.9493), (32.97145850, 0.0076), (33.96786683, 0.0429), (35.96708088, 0.0002)],
             'K': [(38.9637069, 0.932581), (39.96399867, 0.000117), (40.96182597, 0.067302)],
             }
