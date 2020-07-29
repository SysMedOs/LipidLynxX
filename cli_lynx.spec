# -*- mode: python ; coding: utf-8 -*-
import os

block_cipher = None
cwd_path = os.getcwd()

hiddenimports_libs = [
    "aiofiles",
    "click",
    "click-spinner",
    "fastapi",
    "jsonschema",
    "jupyter",
    "ipython",
    "ipywidgets",
    "natsort",
    "nbformat",
    "numpy",
    "openpyxl",
    "pandas",
    "pydantic",
    "pkg_resources",
    "pywin32",
    "python-multipart",
    "regex",
    "starlette",
    "typer",
    "uvicorn",
    "xlrd",
    "xlwt",
]
copy_files = [("lynx/config.ini", "lynx"), ("lynx/static", "lynx/static"), ("lynx/templates", "lynx/templates")]

a = Analysis(
    ["cli_lynx.py"],
    pathex=[cwd_path],
    binaries=[],
    datas=copy_files,
    hiddenimports=hiddenimports_libs,
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="cli_lynx",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    icon="LipidLynxX_Logo.ico",
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="LipidLynxX",
)
