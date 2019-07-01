from dataclasses import dataclass, asdict, field

from epilion.libLION.AbbrElemCalc import ElemCalc

class ElemDict(dict):

    __str_val_keys = ('formula', 'charge_mode', 'adduct')
    __int_val_keys = ('C', 'H', 'O', 'N', 'P', 'S', 'Na', 'K')
    __allowed_keys = ('formula', 'charge_mode', 'adduct', 'C', 'H', 'O', 'N', 'P', 'S', 'Na', 'K')

    __charge_mode_vals = ('positive', 'negative')
    __adducts_vals = {'[M+H]+': 'positive', '[M+NH4]+': 'positive',
                      '[M-H]-': 'negative', '[M+HCOO]-': 'negative', '[M+CH3COO]-': 'negative'}

    def __init__(self, *arg, **kwargs):
        super(ElemDict, self).__init__(*arg, **kwargs)

    def __setitem__(self, key, val):
        if isinstance(key, str):
            if key in self.__allowed_keys:
                if key == 'formula':
                    if not isinstance(val, str):
                        raise ValueError(f'{key} value be an str')
                elif key == 'charge_mode':
                    if not isinstance(val, str):
                        raise ValueError(f'{key} value be an str')
                else:
                    if not isinstance(val, int):
                        raise ValueError(f'{key} value be an int')
            else:
                raise ValueError(f'key must be in list {self.__allowed_keys}')
        else:
            raise ValueError(f'key must be a str in list {self.__allowed_keys}')
        dict.__setitem__(self, key, val)

    def __getattr__(self, key):

        if key in self.__allowed_keys:
            try:
                val = self.__getitem__(key)
            except KeyError:
                if key == 'formula':
                    val = ''
                else:
                    val = 0

            return val
        else:
            raise AttributeError(f'key must be in list {self.__allowed_keys}')

    def __setattr__(self, key, val):
        if key in self.__allowed_keys:
            return self.__setitem__(key, val)
        else:
            raise AttributeError(f'key must be in list {self.__allowed_keys}')


def calc_mz(abbr):
    exactmass = 0.0
    if abbr:
        elem_calc = ElemCalc()
        formula_dct = elem_calc.get_formula(abbr)[1]
        exactmass = elem_calc.get_exactmass(formula_dct)

    return exactmass


@dataclass(frozen=True)
class FormulaDict:

    formula: str = None
    mz: float = field(default_factory=calc_mz(formula))


if __name__ == '__main__':

    fd = ElemDict()

    ffd = FormulaDict(formula='CHO')

    print(asdict(ffd))
    print(ffd.mz)
    ffd.mz = 2.5
    print(ffd.Na)

    fd['formula'] = 'CHO'
    fd['C'] = 10

    print(fd.C)
    fd.C = 1

    print(fd.O)
    # print(fd)
    # fd.bug = 'bug1'
    print(fd)
    fd['BUG'] = 'bug'
    print(fd)
