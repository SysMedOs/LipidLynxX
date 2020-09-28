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

import zmq


def init_client(
    token: str, job: str, data: str, data_type: str = "list", port: int = 2409,
):
    context = zmq.Context()
    # Static analysis of Pycharm/PyLint cannot find runtime-defined names e.g. REQ
    # See: https://github.com/zeromq/pyzmq/issues/1018
    socket = context.socket(zmq.REQ)
    socket.connect(f"tcp://localhost:{port}")
    msg = {
        "token": token,
        "job": job,
        "data": data,
        "data_type": data_type,
    }
    send_msg = json.dumps(msg)
    socket.send(send_msg.encode())
    print(send_msg)
    message = socket.recv()
    print(f"Received reply UUID {token} -  job: {job}, Message: [{message}]")
    socket.close()


def converter_client(token: str, data: str, data_type: str = "list", port: int = 2409):
    init_client(
        token=token, job="convert", data=data, data_type=data_type, port=port,
    )


def equalizer_client(token: str, data: str, data_type: str = "dict", port: int = 2409):
    init_client(
        token=token, job="equalize", data=data, data_type=data_type, port=port,
    )


def linker_client(token: str, data: str, data_type: str = "list", port: int = 2409):
    init_client(
        token=token, job="link", data=data, data_type=data_type, port=port,
    )
