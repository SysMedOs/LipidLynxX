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

from dataclasses import dataclass, asdict, field

from lynx.liblynx.AbbrElemCalc import ElemCalc


class ElemDict(dict):

    __str_val_keys = ("formula", "charge_mode", "adduct")
    __int_val_keys = ("C", "H", "O", "N", "P", "S", "Na", "K")
    __allowed_keys = (
        "formula",
        "charge_mode",
        "adduct",
        "C",
        "H",
        "O",
        "N",
        "P",
        "S",
        "Na",
        "K",
    )

    __charge_mode_vals = ("positive", "negative")
    __adducts_vals = {
        "[M+H]+": "positive",
        "[M+NH4]+": "positive",
        "[M-H]-": "negative",
        "[M+HCOO]-": "negative",
        "[M+CH3COO]-": "negative",
    }

    def __init__(self, *arg, **kwargs):
        super(ElemDict, self).__init__(*arg, **kwargs)

    def __setitem__(self, key, val):
        if isinstance(key, str):
            if key in self.__allowed_keys:
                if key == "formula":
                    if not isinstance(val, str):
                        raise ValueError(f"{key} value must be a str")
                elif key == "charge_mode":
                    if not isinstance(val, str):
                        raise ValueError(f"{key} value must be a str")
                else:
                    if not isinstance(val, int):
                        raise ValueError(f"{key} value must be an int")
            else:
                raise ValueError(f"key must be in list {self.__allowed_keys}")
        else:
            raise ValueError(f"key must be a str in list {self.__allowed_keys}")
        dict.__setitem__(self, key, val)

    def __getattr__(self, key):

        if key in self.__allowed_keys:
            try:
                val = self.__getitem__(key)
            except KeyError:
                if key == "formula":
                    val = ""
                else:
                    val = 0

            return val
        else:
            raise AttributeError(f"key must be in list {self.__allowed_keys}")

    def __setattr__(self, key, val):
        if key in self.__allowed_keys:
            return self.__setitem__(key, val)
        else:
            raise AttributeError(f"key must be in list {self.__allowed_keys}")


def calc_mz(formula):
    exact_mass = 0.0
    print(formula)
    if formula:
        elem_calc = ElemCalc()
        formula_dct = elem_calc.formula_to_elem(formula)
        exact_mass = elem_calc.get_exact_mass(formula_dct)
        print("exact_mass ", exact_mass)

    return exact_mass


@dataclass
class FormulaDict:

    formula: str = None
    exact_mass: float = 0.0

    def __post_init__(self):
        self.__setattr__("exact_mass", calc_mz(self.formula))


if __name__ == "__main__":

    elem_dct = ElemDict()

    usr_formula_dct = FormulaDict(formula="C1H2O3")

    print(asdict(usr_formula_dct))
    print(usr_formula_dct.exact_mass)
    try:
        usr_formula_dct.exact_mass = 2.5
        print(usr_formula_dct.exact_mass)
    except Exception as e:
        print(e)

    print(usr_formula_dct.exact_mass)

    try:
        usr_formula_dct.Na = 1
        print("Na", usr_formula_dct.Na)
    except Exception as e:
        print(e)

    elem_dct["formula"] = "CHO"
    elem_dct["C"] = 10

    print("C", elem_dct.C)
    elem_dct.C = 1

    print("O", elem_dct.O)
    # print(fd)
    # fd.bug = 'bug1'
    print(elem_dct)
    elem_dct["BUG"] = "bug"
    print(elem_dct)
