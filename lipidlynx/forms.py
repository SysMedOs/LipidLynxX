# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
from flask_wtf import FlaskForm
from wtforms import StringField, TextAreaField, FileField
from wtforms.validators import DataRequired, Length


class ConverterTextInputForm(FlaskForm):

    input_id_str = TextAreaField(
        "Paste lipid abbreviations here:", validators=[DataRequired()]
    )


class ConverterTableInputForm(FlaskForm):

    input_file_str = FileField("Select input file:", validators=[DataRequired()])


class EqualizerTableInputForm(FlaskForm):
    input_id_str = TextAreaField(
        "Input LipidLynx level here:", validators=[DataRequired()]
    )
    input_file_str = FileField("Select input file:", validators=[DataRequired()])


class ParserInputForm(FlaskForm):

    lion_id_str = StringField(
        "Paste one lipid abbreviation here:",
        validators=[DataRequired(), Length(max=255)],
    )
