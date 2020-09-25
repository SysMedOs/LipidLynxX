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

from multiprocessing.context import Process

from lynx.models.defaults import default_zmq_worker_runner
from lynx.mq.broker import default_broker
from lynx.mq.worker import general_worker
from lynx.utils.cfg_reader import app_cfg_info
from lynx.utils.ports import check_port


def daemon_lynx(zmq_client_port: int = 2409, zmq_worker_port: int = 2410):
    print("Start ZMQ Broker...")
    Process(target=default_broker, args=(zmq_client_port, zmq_worker_port)).start()
    for w in range(1, default_zmq_worker_runner + 1):
        print(f"Start LipidLynxX ZMQ Worker#{w}...")
        Process(target=general_worker, args=(w, zmq_worker_port)).start()


if __name__ == "__main__":
    # check ports
    checked_zmq_client_port = check_port(
        int(app_cfg_info.get("zmq_client_port", 2409)), task_name="ZMQ client"
    )
    checked_zmq_worker_port = check_port(
        int(app_cfg_info.get("zmq_worker_port", 2410)), task_name="ZMQ worker"
    )
    # run zmq daemon
    daemon_lynx(checked_zmq_client_port, checked_zmq_worker_port)
