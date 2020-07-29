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
import re
from typing import Union, List

import aiohttp
import requests
import natsort
import pandas as pd
import urllib.parse

from lynx.controllers.converter import convert_lipid
from lynx.models.api_models import StyleType
from lynx.models.defaults import kegg_ids, lion_ids

DEFAULT_DB_INFO = {
    "chebi": "https://www.ebi.ac.uk/chebi",
    "hmdb": "https://hmdb.ca",
    "kegg": "https://www.kegg.jp",
    "lion": "http://www.lipidontology.com",
    "lipidbank": "http://lipidbank.jp",
    "lipidmaps": "https://www.lipidmaps.org",
    "pubchem": "https://pubchem.ncbi.nlm.nih.gov",
    "recon2.01": "https://www.vmh.life/#reconmap2",
    "recon3": "https://www.vmh.life/#reconmap",
    "rhea": "https://www.rhea-db.org",
    "swisslipids": "https://www.swisslipids.org",
    "uniportkb": "https://www.uniprot.org",
    "vmh_metabolites": "https://www.vmh.life/#metabolite",
}

CROSS_LINK_DBS = {
    "swisslipids": {
        "chebi": "ChEBI",
        "hmdb": "HMDB",
        "lipidmaps": "LipidMaps",
        "swisslipids": "SwissLipids",
        "rhea": "Rhea",
        "uniportkb": "UniProtKB",
    },
    "lipidmaps": {
        "chebi": "chebi_id",
        "hmdb": "hmdb_id",
        "kegg": "kegg_id",
        "lipidbank": "lipidbank_id",
        "pubchem": "pubchem_cid",
    },
}
CROSS_LINK_APIS = {
    "name_search": {
        "swisslipids": r"https://www.swisslipids.org/api/search?term=<lipid_id>"
    },
    "subclass_search": {
        "lipidmaps": r"https://www.lipidmaps.org/rest/compound/lm_id/<lipid_id>/sub_class",
    },
    "cross_search": {
        "swisslipids": r"https://www.swisslipids.org/api/mapping?from=swisslipids&to=<DB_NAME>&ids=<lipid_id>",
    },
    "link": {
        "chebi": r"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:<lipid_id>",
        "hmdb": r"https://hmdb.ca/metabolites/<lipid_id>",
        "kegg": r"https://www.kegg.jp/dbget-bin/www_bget?cpd:<lipid_id>",
        "lion": r"http://bioportal.bioontology.org/ontologies/LION/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FLION_<lipid_id>",
        "lipidbank": r"http://lipidbank.jp/cgi-bin/detail.cgi?id=<lipid_id>",
        "lipidmaps": r"https://www.lipidmaps.org/data/LMSDRecord.php?LMID=<lipid_id>",
        "pubchem": r"https://pubchem.ncbi.nlm.nih.gov/compound/<lipid_id>",
        "recon2.01": "https://www.vmh.life/minerva/index.xhtml?id=ReconMap-2.01&search=<lipid_id>",
        "recon3": "https://www.vmh.life/minerva/index.xhtml?id=ReconMap-3&search=<lipid_id>",
        "rhea": r"https://www.rhea-db.org/reaction?id=<lipid_id>",
        "swisslipids": r"https://www.swisslipids.org/#/entity/<lipid_id>/",
        "uniportkb": r"https://www.uniprot.org/uniprot/?query=<lipid_id>",
        "vmh_metabolites": "https://www.vmh.life/_api/metabolites/?hmdb=<lipid_id>&format=json",
    },
}

DB_SECTIONS = {
    "Lipid database": ["lipidmaps", "swisslipids"],
    "Metabolites database": ["hmdb", "vmh_metabolites"],
    "Lipid ontology": ["lion"],
    "General database": ["chebi", "pubchem"],
    "Pathways": ["kegg", "recon2.01", "recon3"],
    "Reactions": ["rhea"],
    "Related database": ["uniportkb"],
}


async def get_lmsd_name(lm_id: str = "LMGP01010594") -> str:
    lmsd_base_url = r"https://www.lipidmaps.org/rest/compound/lm_id/<lipid_id>/all"
    ref_url = re.sub(r"<lipid_id>", lm_id, lmsd_base_url)
    lmsd_name = ""
    async with aiohttp.request("GET", ref_url) as r_cross_ref_obj:
        r_cross_ref_status = r_cross_ref_obj.status
        if r_cross_ref_status == 200:
            r_cross_ref_js = await r_cross_ref_obj.json(content_type="application/json")
            lmsd_name = r_cross_ref_js.get("name", "")

    return lmsd_name


async def get_swiss_name(swiss_id: str = "SLM:000000792") -> str:
    swiss_base_url = r"https://www.swisslipids.org/api/entity/<lipid_id>"
    ref_url = re.sub(r"<lipid_id>", swiss_id, swiss_base_url)
    swiss_name = ""
    async with aiohttp.request("GET", ref_url) as r_cross_ref_obj:
        r_cross_ref_status = r_cross_ref_obj.status
        if r_cross_ref_status == 200:
            r_cross_ref_js = await r_cross_ref_obj.json(
                content_type="text/html;charset=utf-8"
            )
            synonyms = r_cross_ref_js.get("synonyms", [])
            for s in synonyms:
                if isinstance(s, dict) and s.get("type", "").lower() == "abbreviation":
                    swiss_name = s.get("name", "")

    return swiss_name


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


def get_lion_id(lipid_name: str = "PC(16:0/20:4)"):
    lion_id = lion_ids.get(lipid_name)
    if lion_id:
        pass
    else:
        lion_id = ""
    return lion_id


def get_lmsd_subclass(lmsd_id: str = "LMGP01010594") -> str:

    sub_class = ""
    if lmsd_id.startswith("LM") and len(lmsd_id) >= 8:
        raw_subclass_str = lmsd_id[2:8]
        if re.match(r"\w\w\d\d\d\d", raw_subclass_str):
            sub_class = raw_subclass_str

    return sub_class


async def get_lmsd_id(
    lipid_name: str = "PC(16:0/18:2(9Z,12Z))", export_url: bool = False
) -> Union[dict, list]:
    lmsd_ids_lst = []
    lmsd_ids_dct = {}
    lmsd_bulk_base_url = "https://www.lipidmaps.org/rest/compound/abbrev/<lipid_id>/all"
    lmsd_chain_base_url = (
        "https://www.lipidmaps.org/rest/compound/abbrev_chains/<lipid_id>/all"
    )
    lipid_name = re.sub(r"/", "_", lipid_name)
    if "_" in lipid_name:
        url_safe_lipid_name = urllib.parse.quote(lipid_name, safe="")
        ref_url = re.sub(r"<lipid_id>", url_safe_lipid_name, lmsd_chain_base_url)
    else:
        url_safe_lipid_name = urllib.parse.quote(lipid_name, safe="")
        ref_url = re.sub(r"<lipid_id>", url_safe_lipid_name, lmsd_bulk_base_url)
    print(ref_url)
    async with aiohttp.request("GET", ref_url) as r_cross_ref_obj:
        r_cross_ref_status = r_cross_ref_obj.status
        if r_cross_ref_status == 200:
            r_cross_ref_js = await r_cross_ref_obj.json(content_type="application/json")
            if isinstance(r_cross_ref_js, dict):
                for lmsd_row in r_cross_ref_js:
                    lmsd_info = r_cross_ref_js.get(lmsd_row)
                    if isinstance(lmsd_info, dict):
                        lmsd_id = lmsd_info.get("lm_id")
                        if lmsd_id:
                            # print(lmsd_id)
                            if export_url:
                                lmsd_ids_dct[lmsd_id] = await get_external_link(
                                    lmsd_id, "lipidmaps"
                                )
                            else:
                                lmsd_ids_lst.append(lmsd_id)
    if export_url:
        lmsd_ids = lmsd_ids_dct
    else:
        lmsd_ids = lmsd_ids_lst
    return lmsd_ids


def get_swiss_id(lipid_name: str = "PC(16:0/18:2(9Z,12Z))") -> requests.Response:
    url_safe_lipid_name = urllib.parse.quote(lipid_name, safe="")
    q_swiss_str = CROSS_LINK_APIS.get("name_search", {}).get("swisslipids")
    q_swiss_str = re.sub(r"<lipid_id>", url_safe_lipid_name, q_swiss_str)
    r_swiss_obj = requests.get(q_swiss_str)

    return r_swiss_obj


async def get_swiss_child(swisslipids_id: str) -> list:
    swiss_children = []
    swiss_base_url = r"https://www.swisslipids.org/api/children?entity_id=<lipid_id>"
    ref_url = re.sub(r"<lipid_id>", str(swisslipids_id), swiss_base_url)
    async with aiohttp.request("GET", ref_url) as r_cross_ref_obj:
        r_cross_ref_status = r_cross_ref_obj.status
        if r_cross_ref_status == 200:
            r_cross_ref_js = await r_cross_ref_obj.json(content_type="text/html")
            if isinstance(r_cross_ref_js, list):
                for swiss_children_info in r_cross_ref_js:
                    if isinstance(swiss_children_info, dict):
                        swiss_children = list(swiss_children_info.keys())

    return swiss_children


async def get_swiss_linked_id(
    swisslipids_id, ref_db, export_url: bool = False
) -> Union[dict, list]:
    ref_dbs = CROSS_LINK_DBS.get("swisslipids", {})
    swiss_base_url = "https://www.swisslipids.org/api/mapping?from=SwissLipids&to=<DB_NAME>&ids=<lipid_id>"
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
            ref_url = re.sub(r"<lipid_id>", str(swisslipids_id), pre_ref_url)
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


async def get_lmsd_linked_ids(
    lm_id: str = "LMGP01010594", ref_dbs: List[str] = None, export_url: bool = False
) -> dict:
    if not ref_dbs:
        ref_dbs = CROSS_LINK_DBS.get("lipidmaps", {})
    lmsd_base_url = r"https://www.lipidmaps.org/rest/compound/lm_id/<lipid_id>/all"
    ref_url = re.sub(r"<lipid_id>", lm_id, lmsd_base_url)
    cross_ref_ids = {}
    async with aiohttp.request("GET", ref_url) as r_cross_ref_obj:
        r_cross_ref_status = r_cross_ref_obj.status
        if r_cross_ref_status == 200:
            r_cross_ref_js = await r_cross_ref_obj.json(content_type="application/json")
            for ref_db in ref_dbs:
                ref_id = r_cross_ref_js.get(ref_dbs.get(ref_db))
                if ref_id:
                    cross_ref_ids[ref_db] = ref_id

    lion_id = get_lion_id(lm_id)
    if lion_id:
        cross_ref_ids["lion"] = lion_id
    if cross_ref_ids:
        pass

    if "kegg" not in cross_ref_ids:
        kegg_id = get_kegg_id(lm_id)
        if kegg_id:
            cross_ref_ids["kegg"] = kegg_id

    return cross_ref_ids


async def get_external_link(ref_id: str, ref_db: str, check_url: bool = False) -> str:

    ref_db_urls = CROSS_LINK_APIS.get("link", {})
    ref_url = ""
    if ref_db in DEFAULT_DB_INFO and ref_id:
        ref_base_url = ref_db_urls.get(ref_db)
        if ref_base_url:
            ref_url = re.sub(r"<lipid_id>", ref_id, ref_base_url)
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
    lipid_name: str = "PC(16:0/18:2(9Z,12Z))",
    export_url: bool = False,
    formatted: bool = True,
) -> dict:
    lipid_name = lipid_name.strip('"')
    linked_ids = {}
    r_swiss_obj = get_swiss_id(lipid_name)
    if r_swiss_obj.status_code == 200:
        r_swiss_js = r_swiss_obj.json()
        for swiss_ref in r_swiss_js:
            swisslipids_name = swiss_ref.get("entity_name")
            swisslipids_id = swiss_ref.get("entity_id")
            if swisslipids_name.startswith("DG"):
                siwss_lv_s_str = re.sub(r"/0:0", "", swisslipids_name)
                siwss_lv_s_str = re.sub(r"\(0:0/", r"(", siwss_lv_s_str)
                siwss_lv_s_str = re.sub(r"/0:0/", r"/", siwss_lv_s_str)
            else:
                siwss_lv_s_str = swisslipids_name
            siwss_lv_s_str = re.sub(
                r"\((\d{1,2}[EZ]?)(,\d{1,2}[EZ]?)+\)", "", siwss_lv_s_str
            )
            siwss_lv_m_str = re.sub(r"/", "_", siwss_lv_s_str)
            print("SWISS_NAMES: ", swisslipids_name, siwss_lv_s_str, siwss_lv_m_str)
            if (
                lipid_name == swisslipids_name
                or lipid_name == siwss_lv_s_str
                or lipid_name == siwss_lv_m_str
            ):
                # linked_ids["swisslipids"] = swisslipids_id
                swiss_ids = [swisslipids_id] + await get_swiss_child(swisslipids_id)
                # print("swiss_ids", swiss_ids)
                for swiss_id in swiss_ids:
                    for cross_ref_db in CROSS_LINK_DBS.get("swisslipids"):
                        cross_ref_ids = await get_swiss_linked_id(
                            swiss_id, cross_ref_db, export_url
                        )
                        if cross_ref_ids:
                            if export_url:
                                existed_links = linked_ids.get(cross_ref_db)
                                if existed_links:
                                    for ref_id in cross_ref_ids:
                                        if ref_id not in existed_links:
                                            existed_links[ref_id] = cross_ref_ids.get(
                                                ref_id
                                            )
                                    linked_ids[cross_ref_db] = existed_links
                                else:
                                    linked_ids[cross_ref_db] = cross_ref_ids
                            else:
                                existed_links = linked_ids.get(cross_ref_db)
                                if existed_links:
                                    existed_links = list(
                                        set(existed_links + cross_ref_ids)
                                    )
                                    linked_ids[cross_ref_db] = existed_links
                                else:
                                    linked_ids[cross_ref_db] = cross_ref_ids
            if linked_ids:
                if "lipidmaps" in linked_ids:
                    lm_ids = linked_ids.get("lipidmaps", [])
                else:
                    lm_ids = await get_lmsd_id(lipid_name, export_url)
                    linked_ids["lipidmaps"] = lm_ids
                if lm_ids:
                    for lm_id in lm_ids:
                        lm_linked_ids = await get_lmsd_linked_ids(lm_id)
                        for ref_db in lm_linked_ids:
                            ref_id = lm_linked_ids.get(ref_db)
                            if export_url:
                                if ref_db in linked_ids:
                                    existed_ids = linked_ids.get(ref_db)
                                    existed_ids[ref_id] = await get_external_link(
                                        ref_id, ref_db
                                    )
                                    linked_ids[ref_db] = existed_ids
                                else:
                                    linked_ids[ref_db] = {
                                        ref_id: await get_external_link(ref_id, ref_db)
                                    }
                            else:
                                if ref_db in linked_ids:
                                    existed_ids = linked_ids.get(ref_db)
                                    existed_ids.append(ref_id)
                                    linked_ids[ref_db] = list(set(existed_ids))
                                else:
                                    linked_ids[ref_db] = [ref_id]
    else:
        pass
    if "hmdb" in linked_ids:
        if export_url:
            hmdb_urls = linked_ids["hmdb"]
            new_hmdb_urls = {}
            for hmdb_id in hmdb_urls:
                if len(hmdb_id) > 9:
                    new_hmdb_urls[hmdb_id] = hmdb_urls[hmdb_id]
            linked_ids["hmdb"] = new_hmdb_urls
        else:
            linked_ids["hmdb"] = [h_id for h_id in linked_ids["hmdb"] if len(h_id) > 9]
    if formatted:
        output_ids = {}
        for ref_group in DB_SECTIONS:
            grouped_dbs = DB_SECTIONS.get(ref_group)
            grouped_ids = {}
            for g_db in grouped_dbs:
                if g_db in linked_ids:
                    grouped_ids[g_db] = linked_ids[g_db]
            output_ids[ref_group] = grouped_ids
        return output_ids
    else:
        return linked_ids


def add_hyperlink(text: str, url: str) -> str:
    return f'=HYPERLINK("{url}", "{text}")'


async def link_lipids(lipid_list: List[str]) -> pd.DataFrame:
    linked_df = pd.DataFrame()
    linked_info_dct = {}
    idx = 1
    for lipid_name in lipid_list:
        if re.match(r"^LM\w\w\d{8}$", lipid_name, re.IGNORECASE):
            safe_lipid_name = await get_lmsd_name(lipid_name)
        elif re.match(r"^SLM:\d{9}$", lipid_name, re.IGNORECASE):
            safe_lipid_name = await get_swiss_name(lipid_name)
        else:
            safe_lipid_name = lipid_name
        search_name = convert_lipid(
            safe_lipid_name, style=StyleType("BracketsShorthand"), level="MAX"
        )
        shorthand_name = convert_lipid(
            safe_lipid_name, style=StyleType("ShorthandNotation"), level="MAX"
        )
        if (
            search_name.startswith("SM(")
            and not search_name.startswith("SM(d")
            and not search_name.startswith("SM(m")
        ):
            search_name = re.sub(r"SM\(", "SM(d", search_name)
        elif (
            search_name.startswith("Cer(")
            and not search_name.startswith("Cer(d")
            and not search_name.startswith("Cer(m")
        ):
            search_name = re.sub(r"Cer\(", "Cer(d", search_name)
        elif (
            search_name.startswith("HexCer(")
            and not search_name.startswith("HexCer(d")
            and not search_name.startswith("HexCer(m")
        ):
            search_name = re.sub(r"HexCer\(", "HexCer(d", search_name)
        elif (
            search_name.startswith("Hex2Cer(")
            and not search_name.startswith("Hex2Cer(d")
            and not search_name.startswith("Hex2Cer(m")
        ):
            search_name = re.sub(r"Hex2Cer\(", "Hex2Cer(d", search_name)
        lynx_name = convert_lipid(
            safe_lipid_name, style=StyleType("LipidLynxX"), level="MAX"
        )
        biopan_name = convert_lipid(safe_lipid_name, style=StyleType("COMP_DB"))
        resources = {
            "Input_name": lipid_name,
            "ShorthandNotation": shorthand_name,
            "LipidLynxX": lynx_name,
            "BioPAN": biopan_name,
        }
        linked_ids = await get_cross_links(
            lipid_name=search_name, export_url=True, formatted=False
        )
        if isinstance(linked_ids, dict):
            for db in linked_ids:
                db_resources = linked_ids.get(db)
                if db_resources and isinstance(db_resources, dict):
                    if len(list(db_resources.keys())) < 2:
                        resources[db] = ";".join(list(db_resources.keys()))
                        resources[f"Link_{db}"] = ";".join(
                            [db_resources.get(i) for i in db_resources]
                        )
                    else:
                        resources[db] = json.dumps(list(db_resources.keys()))
                        resources[f"Link_{db}"] = json.dumps(
                            [db_resources.get(i) for i in db_resources]
                        )
                else:
                    resources[db] = ""
        linked_info_dct[idx] = resources
        print(resources)

        idx += 1
    default_col = ["Input_name", "ShorthandNotation", "LipidLynxX", "BioPAN"]
    if linked_info_dct:
        sum_df = pd.DataFrame(data=linked_info_dct).T
        sum_df_columns = sum_df.columns.tolist()
        link_cols = []
        for col in sum_df_columns:
            if col.startswith("Link_"):
                sum_df_columns.remove(col)
                link_cols.append(col)
            elif col in default_col:
                sum_df_columns.remove(col)
        sum_df_columns = (
            default_col
            + natsort.natsorted(sum_df_columns)
            + natsort.natsorted(link_cols)
        )
        linked_df = pd.DataFrame(sum_df, columns=sum_df_columns)

    return linked_df


if __name__ == "__main__":
    import asyncio

    # # t_id = asyncio.run(get_lmsd_name("LMGP01010594"))
    # t_id = asyncio.run(get_swiss_name("SLM:000000792"))
    # print(t_id)
    # ids = asyncio.run(get_cross_links(t_id, export_url=True, formatted=True))
    # print(ids)
    part = 11
    df = pd.read_excel(f"D:/SysMedOs/LipidLynxX/temp/test_list{part}.xlsx")
    lipid_id_lst = df["ID"].tolist()
    # lipid_id_lst = ["LPC(16:0)", "DG(18:2/20:5)"]
    results_df = asyncio.run(link_lipids(lipid_id_lst))
    results_df.to_excel(f"D:/SysMedOs/LipidLynxX/temp/test_list_part{part}.xlsx")
    print("FIN")
