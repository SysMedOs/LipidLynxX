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

from lynx.models.defaults import default_zmq_client_port


def default_client(data: str, token: str, task: str):
    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.connect(f"tcp://localhost:{default_zmq_client_port}")
    msg = {
        "token": token,
        "job": task,
        "data": data,
    }
    send_msg = json.dumps(msg)
    socket.send(send_msg.encode())
    print(send_msg)
    message = socket.recv()
    print(f"Received reply UUID {token} Message: [{message}]")


def converter_client(token: str, data: str):
    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.connect(f"tcp://localhost:{default_zmq_client_port}")
    msg = {
        "token": token,
        "job": "convert",
        "data": data,
    }
    send_msg = json.dumps(msg)
    socket.send(send_msg.encode())
    message = socket.recv()
    print(f"Received reply UUID {token} Message: [{message}]")
    socket.close()
