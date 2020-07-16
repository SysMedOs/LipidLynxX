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
import requests
import natsort
import urllib.parse

from lynx.models.defaults import kegg_ids

DEFAULT_DB_INFO = {
    "chbei": "https://www.ebi.ac.uk/chebi",
    "hmdb": "https://hmdb.ca",
    "kegg": "https://www.kegg.jp",
    "lion": "http://www.lipidontology.com",
    "lipidbank": "http://lipidbank.jp",
    "lipidmaps": "https://www.lipidmaps.org",
    "pubchem": "https://pubchem.ncbi.nlm.nih.gov",
    "rhea": "https://www.rhea-db.org",
    "swisslipids": "https://www.swisslipids.org",
    "uniportkb": "https://www.uniprot.org",
}

CROSS_LINK_DBS = {
    "swisslipids": {
        "chbei": "ChEBI",
        "hmdb": "HMDB",
        "lipidmaps": "LipidMaps",
        "swisslipids": "SwissLipids",
        "rhea": "Rhea",
        "uniportkb": "UniProtKB",
    },
    "lipidmaps": {
        "chbei": "chebi_id",
        "hmdb": "hmdb_id",
        "kegg": "kegg_id",
        "lipidbank": "lipidbank_id",
        "pubchem": "pubchem_cid",
    },
}
CROSS_LINK_APIS = {
    "name_search": {
        "swisslipids": r"https://www.swisslipids.org/api/search?term=<LIPID_ID>"
    },
    "subclass_search": {
        "lipidmaps": r"https://www.lipidmaps.org/rest/compound/lm_id/<LIPID_ID>/sub_class",
    },
    "cross_search": {
        "swisslipids": r"https://www.swisslipids.org/api/mapping?from=swisslipids&to=<DB_NAME>&ids=<LIPID_ID>",
    },
    "link": {
        "chebi": r"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:<LIPID_ID>",
        "hmdb": r"https://hmdb.ca/metabolites/<LIPID_ID>",
        "kegg": r"https://www.kegg.jp/dbget-bin/www_bget?cpd:<LIPID_ID>",
        "lipidbank": r"http://lipidbank.jp/cgi-bin/detail.cgi?id=<LIPID_ID>",
        "lipidmaps": r"https://www.lipidmaps.org/data/LMSDRecord.php?LMID=<LIPID_ID>",
        "pubchem": r"https://pubchem.ncbi.nlm.nih.gov/compound/<LIPID_ID>",
        "rhea": r"https://www.rhea-db.org/reaction?id=<LIPID_ID>",
        "swisslipids": r"https://www.swisslipids.org/#/entity/<LIPID_ID>/",
        "uniportkb": r"https://www.uniprot.org/uniprot/?query=<LIPID_ID>",
    },
}


def get_swiss_id(lipid_name: str = "PC(16:0/18:2(9Z,12Z))") -> requests.Response:
    url_safe_lipid_name = urllib.parse.quote(lipid_name, safe="")
    q_swiss_str = CROSS_LINK_APIS.get("name_search", {}).get("swisslipids")
    q_swiss_str = re.sub(r"<LIPID_ID>", url_safe_lipid_name, q_swiss_str)
    r_swiss_obj = requests.get(q_swiss_str)

    return r_swiss_obj


def get_lmsd_subclass(lmsd_id: str = "LMGP01010594") -> str:

    sub_class = ""
    if lmsd_id.startswith("LM") and len(lmsd_id) >= 8:
        raw_subclass_str = lmsd_id[2:8]
        if re.match(r"\w\w\d\d\d\d", raw_subclass_str):
            sub_class = raw_subclass_str

    return sub_class


def get_kegg_id(lmsd_id: str = "LMGP01010594"):
    kegg_id = kegg_ids.get(lmsd_id)
    if kegg_id:
        pass
    else:
        lmsd_subclass_id = f"{lmsd_id[:-4]}0000"
        kegg_id = kegg_ids.get(lmsd_subclass_id)
        if kegg_id:
            pass
        else:
            kegg_id = ""
    return kegg_id


async def get_swiss_links(
    swisslipids_id, ref_db, export_url: bool = False
) -> Union[dict, list]:
    ref_dbs = CROSS_LINK_DBS.get("swisslipids", {})
    swiss_base_url = "https://www.swisslipids.org/api/mapping?from=SwissLipids&to=<DB_NAME>&ids=<LIPID_ID>"
    cross_ref_id_urls = {}
    cross_ref_ids = []
    if swisslipids_id and ref_db in ref_dbs:
        if ref_db.lower() == "swisslipids":
            if export_url and swisslipids_id not in cross_ref_id_urls:
                cross_ref_id_urls[swisslipids_id] = await get_external_link(
                    swisslipids_id, ref_db
                )
            else:
                cross_ref_ids.append(swisslipids_id)
        else:
            pre_ref_url = re.sub(r"<DB_NAME>", ref_dbs.get(ref_db), swiss_base_url)
            ref_url = re.sub(r"<LIPID_ID>", str(swisslipids_id), pre_ref_url)
            async with aiohttp.request("GET", ref_url) as r_cross_ref_obj:
                r_cross_ref_status = r_cross_ref_obj.status
                if r_cross_ref_status == 200:
                    r_cross_ref_js = await r_cross_ref_obj.json(
                        content_type="text/html"
                    )
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
                                    if (
                                        export_url
                                        and cross_ref_id not in cross_ref_id_urls
                                    ):
                                        cross_ref_id_urls[
                                            cross_ref_id
                                        ] = await get_external_link(
                                            cross_ref_id, ref_db
                                        )
                                    else:
                                        cross_ref_ids.append(cross_ref_id)
    if export_url:
        cross_ref_id_urls = {
            k: cross_ref_id_urls[k] for k in natsort.natsorted(cross_ref_id_urls.keys())
        }
        return cross_ref_id_urls
    else:
        return list(set(cross_ref_ids))


async def get_external_link(ref_id: str, ref_db: str, check_url: bool = False) -> str:

    ref_db_urls = CROSS_LINK_APIS.get("link", {})
    ref_url = ""
    if ref_id and ref_db in CROSS_LINK_DBS:
        ref_base_url = ref_db_urls.get(ref_db)
        if ref_base_url:
            ref_url = re.sub(r"<LIPID_ID>", ref_id, ref_base_url)
            if check_url:
                print(ref_url)
                async with aiohttp.request("GET", ref_url) as r_cross_ref_obj:
                    r_cross_ref_status = r_cross_ref_obj.status
                    if r_cross_ref_status == 200:
                        pass
                    else:
                        ref_url = ""

    return ref_url


async def get_cross_links(
    lipid_name: str = "PC(16:0/18:2(9Z,12Z))", export_url: bool = False
) -> dict:
    lipid_name = lipid_name.strip('"')
    linked_ids = {}
    r_swiss_obj = get_swiss_id(lipid_name)
    if r_swiss_obj.status_code == 200:
        r_swiss_js = r_swiss_obj.json()
        for swiss_ref in r_swiss_js:
            swisslipids_name = swiss_ref.get("entity_name")
            swisslipids_id = swiss_ref.get("entity_id")
            if swisslipids_name == lipid_name:
                # linked_ids["swisslipids"] = swisslipids_id
                for cross_ref_db in CROSS_LINK_DBS:
                    cross_ref_ids = await get_swiss_links(
                        swisslipids_id, cross_ref_db, export_url
                    )
                    if cross_ref_ids:
                        linked_ids[cross_ref_db] = cross_ref_ids

    else:
        pass

    return linked_ids
