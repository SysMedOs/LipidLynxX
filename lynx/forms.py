# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# LipidLynxX is Dual-licensed
#   For academic and non-commercial use: GPLv2 License:
#   For commercial use: please contact the SysMedOs team by email.
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: lipid annotations converter for large scale lipidomics and epilipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from flask_wtf import FlaskForm
from wtforms import StringField, TextAreaField, FileField
from wtforms.validators import DataRequired, Length


class ConverterTextInputForm(FlaskForm):

    input_id_str = TextAreaField(
        "Paste lipid abbreviations below:", validators=[DataRequired()]
    )


class ConverterTableInputForm(FlaskForm):

    input_file_str = FileField("Select input file:", validators=[DataRequired()])


class ConverterForm(FlaskForm):
    input_id_str = TextAreaField(
        "Paste lipid abbreviations here:", validators=[DataRequired()]
    )
    input_file_str = FileField("Select input file:", validators=[DataRequired()])


class EqualizerInputForm(FlaskForm):
    input_id_str = TextAreaField(
        "Input LipidLynx level(s) here:", validators=[DataRequired()]
    )
    input_file_str = FileField("Select input file:", validators=[DataRequired()])


class ParserInputForm(FlaskForm):

    lion_id_str = StringField(
        "Paste one lipid abbreviation here:",
        validators=[DataRequired(), Length(max=255)],
    )
