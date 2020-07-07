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
from typing import Union

import aiohttp
import urllib.parse
import requests


async def get_cross_ref(lipid_name: str = "PC(16:0/18:2(9Z,12Z))", export_url: bool = False) -> dict:
    lipid_name = lipid_name.strip('"')
    url_safe_lipid_name = urllib.parse.quote(lipid_name, safe="")
    q_swiss_str = f"https://www.swisslipids.org/api/search?term={url_safe_lipid_name}"
    r_swiss_obj = requests.get(q_swiss_str)
    linked_ids = {}
    swiss_cross_refs = ["ChEBI", "HMDB", "LipidMaps", "Rhea", "UniProtKB"]
    if r_swiss_obj.status_code == 200:
        r_swiss_js = r_swiss_obj.json()
        for swiss_ref in r_swiss_js:
            swisslipids_name = swiss_ref.get("entity_name")
            swisslipids_id = swiss_ref.get("entity_id")
            if swisslipids_name == lipid_name:
                linked_ids["SwissLipids"] = swisslipids_id
                for cross_ref_db in swiss_cross_refs:
                    cross_ref_ids = await get_swiss_ref(swisslipids_id, cross_ref_db, export_url)
                    if cross_ref_ids:
                        linked_ids[cross_ref_db] = cross_ref_ids

    else:
        pass

    return linked_ids


async def get_swiss_ref(swisslipids_id, ref_db, export_url: bool = False) -> Union[dict, list]:

    ref_dbs = ["ChEBI", "HMDB", "LipidMaps", "Rhea", "UniProtKB"]

    swiss_base_url = "https://www.swisslipids.org/api/mapping?from=SwissLipids&to=REF_DB_NAME&ids=LIPID_ID"
    cross_ref_id_urls = {}
    cross_ref_ids = []
    if swisslipids_id and ref_db in ref_dbs:
        pre_ref_url = re.sub(r"REF_DB_NAME", str(ref_db), swiss_base_url)
        ref_url = re.sub(r"LIPID_ID", str(swisslipids_id), pre_ref_url)
        async with aiohttp.request('GET', ref_url) as r_cross_ref_obj:
            r_cross_ref_status = r_cross_ref_obj.status
            if r_cross_ref_status == 200:
                r_cross_ref_js = await r_cross_ref_obj.json(content_type='text/html')
                for cross_ref_info in r_cross_ref_js:
                    swiss_id = cross_ref_info.get("from", [{}]).get("id")
                    cross_ref_id_lst = cross_ref_info.get("to")
                    if (
                            swiss_id == swisslipids_id
                            and isinstance(cross_ref_id_lst, list)
                            and cross_ref_id_lst
                    ):
                        for cross_ref_id_info in cross_ref_id_lst:
                            cross_ref_id = str(cross_ref_id_info.get("id"))
                            if cross_ref_id:
                                if export_url and cross_ref_id not in cross_ref_id_urls:
                                    cross_ref_id_urls[cross_ref_id] = await get_ref_link(
                                        cross_ref_id, ref_db
                                    )
                                else:
                                    cross_ref_ids.append(cross_ref_id)
    if export_url:
        return cross_ref_id_urls
    else:
        return list(set(cross_ref_ids))


async def get_ref_link(ref_id: str, ref_db: str) -> str:

    ref_dbs = ["ChEBI", "HMDB", "LipidMaps", "Rhea", "SwissLipids", "UniProtKB"]

    ref_db_urls = {
        "ChEBI": r"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:LIPID_ID",
        "HMDB": r"https://hmdb.ca/metabolites/LIPID_ID",
        "LipidMaps": r"https://www.lipidmaps.org/data/LMSDRecord.php?LMID=LIPID_ID",
        "Rhea": r"https://www.rhea-db.org/reaction?id=LIPID_ID",
        "UniProtKB": r"https://www.uniprot.org/uniprot/?query=LIPID_ID",
    }
    ref_url = ""
    if ref_id and ref_db in ref_dbs:
        ref_base_url = ref_db_urls.get(ref_db)
        if ref_base_url:
            ref_url = re.sub(r"LIPID_ID", ref_id, ref_base_url)
            # print(ref_url)
            # async with aiohttp.request('GET', ref_url) as r_cross_ref_obj:
            #     r_cross_ref_status = r_cross_ref_obj.status
            #     if r_cross_ref_status == 200:
            #         valid_ref_url = ref_url

    return ref_url
