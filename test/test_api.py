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

import json
import os
import sys
import unittest

lynx_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, lynx_path + "/../")

from lynx import app
from lynx.utils.log import app_logger


# r = requests.get(api_url, params={"data": json.dumps(data)}).json()

#
# class APITest(unittest.TestCase):
#     def setUp(self):
#         logger.debug("SETUP APITest...")
#         self.app = app.test_client()
#
#     def test_convert(self):
#         api_url = api_url_info.get("convert")
#         logger.info(f"Test API: {api_url}")
#         test_data = [
#             "GM3(d18:1/18:2(9Z,11Z)(12OH))",
#             "TG P-18:1_18:2(9Z,11Z)(12OH)_18:1(9)(11OH)",
#             "CL(1'-[18:1(9Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",
#             "TG(16:0/18:2/9:0<oxo{9}>)",
#         ]
#         data_json = {"data": json.dumps(test_data)}
#         logger.info(f"Data for API: {data_json}")
#         response = self.app.get(
#             api_url, headers={"Content-Type": "application/json"}, data=data_json
#         )
#         if response:
#             logger.info(f"Test output: {response}")
#             logger.info(f"Test passed: {api_url}")
#
#     def test_convert_string(self):
#         api_url = api_url_info.get("convert_str")
#         logger.info(f"Test API: {api_url}")
#         test_data = "TG P-18:1_18:2(9Z,11Z)(12OH)_18:1(9)(11OH)"
#         data_json = {"data": json.dumps(test_data)}
#         logger.info(f"Data for API: {data_json}")
#         response = self.app.get(
#             api_url, headers={"Content-Type": "application/json"}, data=data_json
#         )
#         if response:
#             logger.info(f"Test output: {response}")
#             logger.info(f"Test passed: {api_url}")
#
#     def test_convert_list(self):
#         api_url = api_url_info.get("convert_list")
#         logger.info(f"Test API: {api_url}")
#         test_data = [
#             "GM3(d18:1/18:2(9Z,11Z)(12OH))",
#             "TG P-18:1_18:2(9Z,11Z)(12OH)_18:1(9)(11OH)",
#             "CL(1'-[18:1(9Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",
#             "TG(16:0/18:2/9:0<oxo{9}>)",
#         ]
#         data_json = {"data": json.dumps(test_data)}
#         logger.info(f"Data for API: {data_json}")
#         response = self.app.get(
#             api_url, headers={"Content-Type": "application/json"}, data=data_json
#         )
#         if response:
#             logger.info(f"Test output: {response}")
#             logger.info(f"Test passed: {api_url}")
#
#     def test_convert_dict(self):
#         api_url = api_url_info.get("convert_dict")
#         logger.info(f"Test API: {api_url}")
#         test_data = {
#             "TestInput": [
#                 "GM3(d18:1/18:2(9Z,11Z)(12OH))",
#                 "TG P-18:1_18:2(9Z,11Z)(12OH)_18:1(9)(11OH)",
#                 "CL(1'-[18:1(9Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",
#                 "TG(16:0/18:2/9:0<oxo{9}>)",
#             ]
#         }
#         data_json = {"data": json.dumps(test_data)}
#         logger.info(f"Data for API: {data_json}")
#         response = self.app.get(
#             api_url, headers={"Content-Type": "application/json"}, data=data_json
#         )
#         if response:
#             logger.info(f"Test output: {response}")
#             logger.info(f"Test passed: {api_url}")
#
#     def tearDown(self):
#         pass
