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


class ElemCalc:
    def __init__(self):

        pa_hg_elem = {"C": 0, "H": 3, "O": 4, "P": 1, "N": 0}
        pc_hg_elem = {"C": 5, "H": 14, "O": 4, "P": 1, "N": 1}
        pe_hg_elem = {"C": 2, "H": 8, "O": 4, "P": 1, "N": 1}
        pg_hg_elem = {"C": 3, "H": 9, "O": 6, "P": 1, "N": 0}
        pi_hg_elem = {"C": 6, "H": 13, "O": 9, "P": 1, "N": 0}
        pip_hg_elem = {"C": 6, "H": 14, "O": 12, "P": 2, "N": 0}
        ps_hg_elem = {"C": 3, "H": 8, "O": 6, "P": 1, "N": 1}
        tg_hg_elem = {"C": 0, "H": 0, "O": 0, "P": 0, "N": 0}
        fa_hg_elem = {"C": 0, "H": 0, "O": 0, "P": 0, "N": 0}

        self.lipid_hg_elem_dct = {
            "PA": pa_hg_elem,
            "PC": pc_hg_elem,
            "PE": pe_hg_elem,
            "PG": pg_hg_elem,
            "PI": pi_hg_elem,
            "PS": ps_hg_elem,
            "PIP": pip_hg_elem,
            "TG": tg_hg_elem,
            "FA": fa_hg_elem,
        }

        self.glycerol_bone_elem_dct = {"C": 3, "H": 2}
        self.link_o_elem_dct = {"O": -1, "H": 2}
        self.link_p_elem_dct = {"O": -1}

        self.periodic_table_dct = {
            "H": [(1.0078250321, 0.999885), (2.0141017780, 0.0001157)],
            "D": [(2.0141017780, 0.0001157)],
            "C": [(12.0, 0.9893), (13.0033548378, 0.0107)],
            "N": [(14.0030740052, 0.99632), (15.0001088984, 0.00368)],
            "O": [
                (15.9949146221, 0.99757),
                (16.99913150, 0.00038),
                (17.9991604, 0.00205),
            ],
            "Na": [(22.98976967, 1.0)],
            "P": [(30.97376151, 1.0)],
            "S": [
                (31.97207069, 0.9493),
                (32.97145850, 0.0076),
                (33.96786683, 0.0429),
                (35.96708088, 0.0002),
            ],
            "K": [
                (38.9637069, 0.932581),
                (39.96399867, 0.000117),
                (40.96182597, 0.067302),
            ],
        }

        self.elem_rgx = re.compile(r"([CDHNOPKS]|Na)(\d{0,3})")
        self.elem_range_check_rgx = re.compile(r"\d{4,}")

    @staticmethod
    def decode_abbr(abbr):

        pl_checker = re.compile(r"(P[ACEGSI])([(])(.*)([)])")
        pip_checker = re.compile(r"(PIP)([(])(.*)([)])")
        tg_checker = re.compile(r"(TG)([(])(.*)([)])")
        fa_checker = re.compile(r"(FA)(\d{1,2})([:])(\d{1,2})")
        fa_short_checker = re.compile(r"(\d{1,2})([:])(\d{1,2})")
        fa_o_checker = re.compile(r"(O-)(\d{1,2})([:])(\d)")
        fa_p_checker = re.compile(r"(P-)(\d{1,2})([:])(\d)")

        # Check PL Type
        _pl_typ = ""
        bulk_fa_typ = ""
        bulk_fa_linker = ""
        bulk_fa_c = 0
        bulk_fa_db = 0
        lyso_fa_linker_dct = {"sn1": "", "sn2": ""}

        if pl_checker.match(abbr):
            # print('PL')
            pl_re_chk = pl_checker.match(abbr)
            pl_typ_lst = pl_re_chk.groups()
            _pl_typ = pl_typ_lst[0]
            bulk_fa_typ = pl_typ_lst[2]
        if pip_checker.match(abbr):
            # print('PIP')
            pip_re_chk = pip_checker.match(abbr)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            bulk_fa_typ = pip_typ_lst[2]
        if tg_checker.match(abbr):
            # print('TG')
            tg_re_chk = tg_checker.match(abbr)
            tg_typ_lst = tg_re_chk.groups()
            _pl_typ = tg_typ_lst[0]
            bulk_fa_typ = tg_typ_lst[2]
        if fa_checker.match(abbr):
            # print('FA')
            _pl_typ = "FA"
            bulk_fa_typ = abbr
            fa_chk = fa_checker.match(abbr)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
            bulk_fa_linker = "A-"
            lyso_fa_linker_dct = {"A": ""}
        if fa_short_checker.match(abbr):
            # print('FA')
            _pl_typ = "FA"
            bulk_fa_typ = abbr
            fa_chk = fa_checker.match(abbr)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[0]
            bulk_fa_db = bulk_fa_lst[2]
            bulk_fa_linker = "A-"
            lyso_fa_linker_dct = {"A": ""}
        if fa_o_checker.match(abbr):
            # print('FA')
            _pl_typ = "FA"
            bulk_fa_typ = abbr
            fa_chk = fa_o_checker.match(abbr)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
            bulk_fa_linker = "O-"
            lyso_fa_linker_dct = {"O": ""}
        if fa_p_checker.match(abbr):
            # print('FA')
            _pl_typ = "FA"
            bulk_fa_typ = abbr
            fa_chk = fa_p_checker.match(abbr)
            bulk_fa_lst = fa_chk.groups()
            bulk_fa_c = bulk_fa_lst[1]
            bulk_fa_db = bulk_fa_lst[3]
            bulk_fa_linker = "P-"
            lyso_fa_linker_dct = {"P": ""}

        if _pl_typ in ["PL", "PA", "PC", "PE", "PG", "PI", "PS"]:
            if fa_short_checker.match(bulk_fa_typ):
                bulk_fa_linker = "A-A-"
                lyso_fa_linker_dct = {"A": ""}
                fa_chk = fa_short_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[0]
                bulk_fa_db = bulk_fa_lst[2]
            # elif fa_short_checker.match(bulk_fa_typ):
            #     bulk_fa_linker = ''
            #     lyso_fa_linker_dct = {'A': ''}
            #     fa_chk = fa_short_checker.match(bulk_fa_typ)
            #     bulk_fa_lst = fa_chk.groups()
            #     bulk_fa_c = bulk_fa_lst[0]
            #     bulk_fa_db = bulk_fa_lst[2]
            elif fa_o_checker.match(bulk_fa_typ):
                bulk_fa_linker = "O-A-"
                lyso_fa_linker_dct = {
                    "O": "",
                    "A": "O-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_o_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]
            elif fa_p_checker.match(bulk_fa_typ):
                bulk_fa_linker = "P-A-"
                lyso_fa_linker_dct = {
                    "P": "",
                    "A": "P-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_p_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]
        elif _pl_typ in ["TG"]:
            if fa_checker.match(bulk_fa_typ):
                bulk_fa_linker = "A-A-A-"
                lyso_fa_linker_dct = {"A": ""}
                fa_chk = fa_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[0]
                bulk_fa_db = bulk_fa_lst[2]
            elif fa_o_checker.match(bulk_fa_typ):
                bulk_fa_linker = "O-A-A-"
                lyso_fa_linker_dct = {
                    "O": "",
                    "A": "O-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_o_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]
            elif fa_p_checker.match(bulk_fa_typ):
                bulk_fa_linker = "P-A-A-"
                lyso_fa_linker_dct = {
                    "P": "",
                    "A": "P-",
                }  # link of the other sn after NL of this sn
                fa_chk = fa_p_checker.match(bulk_fa_typ)
                bulk_fa_lst = fa_chk.groups()
                bulk_fa_c = bulk_fa_lst[1]
                bulk_fa_db = bulk_fa_lst[3]

        bulk_fa_c = int(bulk_fa_c)
        bulk_fa_db = int(bulk_fa_db)

        lipid_info_dct = {
            "TYPE": _pl_typ,
            "LINK": bulk_fa_linker,
            "C": bulk_fa_c,
            "DB": bulk_fa_db,
            "LYSO_LINK": lyso_fa_linker_dct,
        }

        return lipid_info_dct

    def get_neutral_elem(self, abbr):

        usr_lipid_info_dct = self.decode_abbr(abbr)

        lipid_type = usr_lipid_info_dct["TYPE"]

        if lipid_type in self.lipid_hg_elem_dct.keys():
            if lipid_type in ["FA"]:
                # print(abbr)
                tmp_lipid_elem_dct = {
                    "C": usr_lipid_info_dct["C"],
                    "O": 2,
                    "H": (usr_lipid_info_dct["C"] * 2 - usr_lipid_info_dct["DB"] * 2),
                }
                if usr_lipid_info_dct["LINK"] == "":
                    pass
                elif usr_lipid_info_dct["LINK"] == "O-":
                    tmp_lipid_elem_dct["O"] += -1
                    tmp_lipid_elem_dct["H"] += 2
                elif usr_lipid_info_dct["LINK"] == "P-":
                    tmp_lipid_elem_dct["O"] += -1
                else:
                    pass

                return tmp_lipid_elem_dct

            if lipid_type in ["PA", "PC", "PE", "PG", "PI", "PIP", "PS"]:
                tmp_lipid_elem_dct = self.lipid_hg_elem_dct[
                    usr_lipid_info_dct["TYPE"]
                ].copy()
                tmp_lipid_elem_dct["O"] += 4
                tmp_lipid_elem_dct["C"] += (
                    self.glycerol_bone_elem_dct["C"] + usr_lipid_info_dct["C"]
                )
                tmp_lipid_elem_dct["H"] += (
                    self.glycerol_bone_elem_dct["H"]
                    + usr_lipid_info_dct["C"] * 2
                    - usr_lipid_info_dct["DB"] * 2
                )  # DBE = DB + 2xC=O from FA
                if usr_lipid_info_dct["LINK"] == "O-A-":
                    tmp_lipid_elem_dct["O"] += -1
                    tmp_lipid_elem_dct["H"] += 2
                elif usr_lipid_info_dct["LINK"] == "P-A-":
                    tmp_lipid_elem_dct["O"] += -1
                else:
                    pass

                return tmp_lipid_elem_dct

            elif lipid_type in ["TG"]:
                tmp_lipid_elem_dct = self.lipid_hg_elem_dct[
                    usr_lipid_info_dct["TYPE"]
                ].copy()
                tmp_lipid_elem_dct["O"] += 6
                tmp_lipid_elem_dct["C"] += (
                    self.glycerol_bone_elem_dct["C"] + usr_lipid_info_dct["C"]
                )
                tmp_lipid_elem_dct["H"] += (
                    self.glycerol_bone_elem_dct["H"]
                    + usr_lipid_info_dct["C"] * 2
                    - usr_lipid_info_dct["DB"] * 2
                )  # DBE = DB + 2xC=O from FA
                if usr_lipid_info_dct["LINK"] == "O-A-A-":
                    tmp_lipid_elem_dct["O"] += -1
                    tmp_lipid_elem_dct["H"] += 2
                elif usr_lipid_info_dct["LINK"] == "P-A-A-":
                    tmp_lipid_elem_dct["O"] += -1

                return tmp_lipid_elem_dct

            else:
                return {"C": 0, "H": 0, "O": 0, "P": 0}
        else:
            return {"C": 0, "H": 0, "O": 0, "P": 0}

    @staticmethod
    def get_elem_with_charge(elem_dct, charge="[M-H]-"):

        chg_elem_dct = elem_dct.copy()
        if charge == "[M-H]-":
            chg_elem_dct["H"] += -1
        elif charge == "[M+HCOO]-" or charge == "[M+FA-H]-":
            chg_elem_dct["H"] += 1
            chg_elem_dct["C"] += 1
            chg_elem_dct["O"] += 2
        elif charge == "[M+CH3COO]-":
            chg_elem_dct["H"] += 3
            chg_elem_dct["C"] += 2
            chg_elem_dct["O"] += 2
        elif charge == "[M+OAc]-":
            chg_elem_dct["H"] += 3
            chg_elem_dct["C"] += 2
            chg_elem_dct["O"] += 2
        elif charge == "[M+H]+":
            chg_elem_dct["H"] += 1
        elif charge == "[M+NH4]+":
            chg_elem_dct["H"] += 4
            if "N" in chg_elem_dct:
                chg_elem_dct["N"] += 1
            else:
                chg_elem_dct["N"] = 1
        elif charge in ["[M+H-H2O]+", "[M-H2O+H]+"]:
            chg_elem_dct["O"] += -1
            chg_elem_dct["H"] += -1
        elif charge == "[M+Na]+":
            if "Na" in chg_elem_dct:
                chg_elem_dct["Na"] += 1
            else:
                chg_elem_dct["Na"] = 1
        elif charge == "[M+K]+":
            if "K" in chg_elem_dct:
                chg_elem_dct["K"] += 1
            else:
                chg_elem_dct["K"] = 1

        return chg_elem_dct

    def get_charged_elem(self, abbr, charge="[M-H]-"):

        neutral_elem_dct = self.get_neutral_elem(abbr)
        lipid_elem_dct = self.get_elem_with_charge(neutral_elem_dct, charge)

        return lipid_elem_dct

    @staticmethod
    def elem_to_formula(elem_dct, charge="[M-H]-"):

        formula_str = "C{c}H{h}".format(c=elem_dct["C"], h=elem_dct["H"])

        for elem in ["N", "O", "P", "S", "Na", "K"]:

            if elem in elem_dct:
                if elem_dct[elem] == 1:
                    formula_str += elem
                elif elem_dct[elem] > 1:
                    formula_str += f"{elem}{elem_dct[elem]}"

        if charge in ["neutral", "Neutral", "", None]:
            pass
        elif charge in ["[M-H]-", "[M+HCOO]-", "[M+CH3HCOO]-"]:
            formula_str += "-"
        elif charge in [
            "[M+H]+",
            "[M+NH4]+",
            "[M+Na]+",
            "[M+K]+",
            "[M+H-H2O]+",
            "[M-H2O+H]+",
        ]:
            formula_str += "+"

        return formula_str

    def formula_to_elem(self, formula):

        if re.search(self.elem_range_check_rgx, formula):
            raise ValueError(
                f"Can not parse formula {formula}: number of each element should < 999"
            )

        f_found_lst = re.findall(self.elem_rgx, formula)
        re_str = "".join(["".join(r) for r in f_found_lst])
        if re_str != formula:
            print(f"{re_str} != {formula}")
            raise ValueError(
                f"Some elements not supported. Input-> {formula} != {re_str} <- Parsed"
            )

        elem_dct = {}
        if f_found_lst:
            for f_match in f_found_lst:
                if f_match[0] in self.periodic_table_dct:
                    if f_match[1]:
                        try:
                            elem_dct[f_match[0]] = int(f_match[1])
                        except ValueError:
                            raise ValueError
                    else:
                        elem_dct[f_match[0]] = 1
                else:
                    raise ValueError(
                        f"Can not parse formula {formula}: Some elements not supported."
                    )
        else:
            raise ValueError

        return elem_dct

    def get_formula(self, abbr, charge=""):

        if charge in ["neutral", "Neutral", "", None]:

            elem_dct = self.get_neutral_elem(abbr)
        else:
            elem_dct = self.get_charged_elem(abbr, charge=charge)

        formula_str = self.elem_to_formula(elem_dct, charge)

        return formula_str, elem_dct

    def get_exact_mass(self, elem_dct):

        mono_mz = 0.0
        for _elem in elem_dct.keys():
            mono_mz += elem_dct[_elem] * self.periodic_table_dct[_elem][0][0]

        return round(mono_mz, 6)


if __name__ == "__main__":

    usr_bulk_abbr_lst = ["TG(P-48:2)", "PC(O-36:3)", "PC(P-36:3)", "PC(36:3)"]
    charge_lst = ["[M+NH4]+", "[M-H]-", "[M+HCOO]-", "[M+OAc]-"]
    # usr_bulk_abbr_lst = ['PC(36:3)', 'PC(O-36:3)', 'PC(P-36:3)']
    # charge_lst = ['', '[M-H]-', '[M+HCOO]-', '[M+OAc]-']

    abbr2formula = ElemCalc()

    for usr_abbr in usr_bulk_abbr_lst:
        for _charge in charge_lst:
            usr_formula, usr_elem_dct = abbr2formula.get_formula(
                usr_abbr, charge=_charge
            )
            print(usr_abbr, _charge)
            print(usr_elem_dct)
            print(usr_formula)

    f_lst = ["C48H93NO11S", "C148H93NO11S", "C1148H93NO11S", "XC148H93NO11S"]

    for f in f_lst:
        try:
            print(f)
            print(abbr2formula.formula_to_elem(f))
        except Exception as e:
            print(e)
