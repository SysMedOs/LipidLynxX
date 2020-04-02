import json
from typing import List

from jsonschema import Draft7Validator

from lynx.utils.log import logger


def seg_to_str(in_list: List[str], sep: str = ",") -> str:
    """
    Combine multiple segments into one str without space and separator in the end
    Args:
        in_list: input list to be joined.
        sep: separator used to join str segments, default value is ","

    Returns:
        out_str: the joined str

    """
    in_list = filter(None, in_list)
    out_str = sep.join(in_list)
    out_str = out_str.strip(sep)
    if out_str is None:
        out_str = ""
    return out_str


def check_json(validator: Draft7Validator, json_obj: json) -> bool:
    is_valid = False
    if validator.is_valid(json_obj):
        logger.debug(f"JSON Schema check PASSED.")
        is_valid = True
    else:
        for e in validator.iter_errors(json_obj):
            logger.error(e)
        raise Exception(f"JSON Schema check FAILED. Json string: {json_obj}")
    return is_valid


def clean_dct(dct: dict) -> dict:

    if dct:
        for k in dct:
            v = dct[k]
            if isinstance(v, str):
                dct[k] = [v]
            elif isinstance(v, list):
                dct[k] = list(filter(None, v))
            else:
                pass
    else:
        pass

    return dct
