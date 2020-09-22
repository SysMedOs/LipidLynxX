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

from lynx.mq.broker import default_broker
from lynx.mq.worker import converter_worker


def daemon_lynx():

    print("Start ZMQ Broker...")
    Process(target=default_broker).start()
    for w in range(1, 5):
        print(f"Start LipidLynxX ZMQ Worker#{w}...")
        Process(target=converter_worker, args=(w,)).start()


if __name__ == "__main__":
    daemon_lynx()
