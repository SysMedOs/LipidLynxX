# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

if __name__ == "__main__":

    import webbrowser
    from lynx import app
    from lynx.models.defaults import cfg_info_dct

    base_url = cfg_info_dct.get("base_url", "http://127.0.0.1:5000")
    webbrowser.open(f"{base_url}/lynx", new=1, autoraise=True)
    app.run(debug=True)
