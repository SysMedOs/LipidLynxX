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

from multiprocessing import Process
import webbrowser
import uvicorn

from lynx.app import app, app_url, app_port
from lynx.tasks.broker import default_broker
from lynx.tasks.worker import default_worker, converter_worker


def start_lynx():
    print("Start Browser: ")
    webbrowser.open(f"http://{app_url}:{app_port}", new=1, autoraise=True)
    print("Start LipidLynxX ZMQ Broker: ")
    Process(target=default_broker).start()
    for w in range(1, 5):
        print(f"Start LipidLynxX ZMQ Worker#{w}: ")
        Process(target=converter_worker, args=(w,)).start()
    print("Start LipidLynxX Main Application: ")
    uvicorn.run(app, host=app_url, port=app_port)
    # print("Start LipidLynxX Main Application: ")
    # Process(
    #     target=uvicorn.run, args=(app,), kwargs={"host": app_url, "port": app_port}
    # ).start()


if __name__ == "__main__":
    start_lynx()
