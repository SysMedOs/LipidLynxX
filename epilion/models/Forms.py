# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
from flask_wtf import Form
from wtforms import StringField, TextAreaField
from wtforms.validators import DataRequired, Length


class ConverterInputForm(Form):

    input_id_str = TextAreaField(
        "Paste lipid abbreviations here:", validators=[DataRequired()]
    )


class ParserInputForm(Form):

    lion_id_str = StringField(
        "Paste one lipid abbreviation here:",
        validators=[DataRequired(), Length(max=255)],
    )
