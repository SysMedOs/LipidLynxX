# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

if __name__ == "__main__":

    import webbrowser
    from lynx import app

    webbrowser.open("http://127.0.0.1:5000/lynx", new=1, autoraise=True)
    app.run()
