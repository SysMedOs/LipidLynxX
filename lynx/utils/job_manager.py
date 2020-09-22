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
import hashlib
import os
import time
from typing import Union

from lynx.models.api_models import JobType
from lynx.models.defaults import default_temp_folder


def create_job_token(job_type: JobType) -> str:
    time_tag = str(time.time())
    # print(job_type.job, type(job_type.job))
    if not isinstance(job_type, JobType):
        raise TypeError(f"Required type: JobType, received type: {type(job_type)}")
    code_to_hash = "@".join([job_type.job, time_tag])
    hash_obj = hashlib.sha1(code_to_hash.encode())
    job_token = hash_obj.hexdigest()

    return job_token


def save_job(job_token: str, data: Union[str, dict]):
    if len(job_token) > 5:
        job_file = os.path.join(default_temp_folder, f"{job_token}.json")
        print(f"Try to save: {job_file}")
        if isinstance(data, dict):
            with open(job_file, mode="w", encoding="utf-8") as s_obj:
                json.dump(data, s_obj)
        else:
            try:
                data_dct = json.loads(data)
                if data_dct:
                    with open(job_file, mode="w", encoding="utf-8") as s_obj:
                        json.dump(data_dct, s_obj)
            except json.decoder.JSONDecodeError:
                raise ValueError("data must be JSON compatible.")
    else:
        raise ValueError("Invalid session id.")


def is_job_finished(job_token: str) -> bool:

    job_file = os.path.join(default_temp_folder, f"{job_token}.json")
    if os.path.isfile(job_file):
        return True
    else:
        return False


def get_job_output(job_token: str) -> dict:
    job_output = {}
    if is_job_finished(job_token):
        job_file = os.path.join(default_temp_folder, f"{job_token}.json")
        with open(job_file, mode="r", encoding="utf-8") as s_obj:
            job_output = json.load(s_obj)

    return job_output


if __name__ == "__main__":
    # test_dct = {"x": 123}
    # test_text = json.dumps(test_dct)
    # test_js = json.loads(test_text)
    # s1_id = create_job_token(test_dct)
    # save_session(s1_id, test_dct)
    # s2_id = create_job_token(test_text)
    # save_session(s2_id, test_text)
    test_dct = JobType(job="convert")
    s1_id = create_job_token(test_dct)
    # save_session(s1_id, test_dct)

    print(s1_id)
    print(get_job_output(s1_id))

    # print(s2_id)
    # print(get_job_output(s2_id))
