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

import socket


def try_port(port: int) -> bool:
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = False
    try:
        sock.bind(("0.0.0.0", port))
        result = True
    except Exception as e:
        print(f"{e}")
        print(f"Port: [{port}] already in use, try to find next available port.")
    sock.close()
    return result


def check_port(port: int, task_name: str = None) -> int:
    is_port_available = False
    tmp_port = port
    if task_name:
        task_info = f" for {task_name}"
    else:
        task_info = ""
    while is_port_available is False:
        is_port_available = try_port(tmp_port)
        if is_port_available:
            if tmp_port == port:
                print(f"Port: [{port}] is available{task_info}.")
            else:
                print(
                    f"Port: [{port}] already in use, switch to next available port: [{tmp_port}]{task_info}."
                )
        else:
            tmp_port += 3

    return tmp_port


if __name__ == "__main__":
    usr_port = 1399
    available_port = check_port(usr_port, "LipidLynxX")
    print(available_port)
