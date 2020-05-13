# LipidLynxX 

![LipidLynx_Logo](doc/images/LipidLynxX_Logo_128.jpg)

![Platforms](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-blue.svg)
![GitHub (Pre-)Release Date](https://img.shields.io/github/release-date-pre/SysMedOs/LipidLynxX.svg)
![total downloads](https://img.shields.io/github/downloads/SysMedOs/LipidLynxX/total.svg?color=orange)
![GitHub last commit](https://img.shields.io/github/last-commit/SysMedOs/LipidLynxX.svg)

The LipidLynxX project is aimed to provide a unified identifier for major lipids, especially oxidized lipids
in the epilipidome.

![LipidLynx_01_Home](doc/images/LipidLynxX_Start_Chromium.png)

## Try LipidLynxX simple converter demo on [`mybinder.org`](https://mybinder.org)  🆕 

**This demo is always updated automatically to the latest source code on the master branch.**
To preview the latest changes on the converter without dealing with source code.

Just click this button 👉
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ZhixuNi/LipidLynxX/master?filepath=converter_notebook.ipynb)

And wait a bit ☕ Binder and Jupyter Notebook will prepare LipidLynxX demo for you.

- You can paste a list of lipid abbreviations and download the output table as `.csv` or `.xlsx`.

- If you observed some IDs not converted in the Windows .exe version, try this demo to see if it got fixed.

- You can run the notebook named `converter_notebook.ipynb` in this repository as well.

## Important Notice

The current LipidLynxX source code was tested using our collection of lipid abbreviations
for major lipid classes from following databases and programs:

- Databases (5):
  - `HMDB`, `LIPID MAPS LMSD & COMP_DB`, `LipidHome`, `RefMet`, `SwissLipids`

- Programs (17):

  - `ALEX123 lipid calculator`, `Greazy`, `LDA 2`, `LipidBlast`, `LipidCreator`, `LipiDex`, `LipidFrag`, `LipidHunter`,
      `LipidMatch`, `LipidPro`, `LipidSearch`, `Lipostar`, `LIQUID`, `LPPtiger`, `MetFrag`, `MS-DIAL`, `MZmine2`

- Common abbreviations (customizable):
  -  Abbreviations such as DHA, PAPE, PLPC, PONPC .etc are also included as `defined alias`.
  detailed settings can be found in `lynx/configurations/defined_alias.json`

**If your database / program is not included in the list above**, you can test if any of the configuration files located in `lynx/configurations/rules/input` would fit to your database / program.
If conversion is not possible, please contact us so that we can help you to generate suitable configuration file.

A robust and accurate converter can only be achieved by community-wide collaborations, thus any issue reports from general users and developers are welcome and will improve LipidLynxX project.

Thus, if you meet any issues during using LipidLynxX, please [report your issue
here](https://github.com/SysMedOs/LipidLynxX/issues)

### Notice to general users

An easy to use .exe version for Windows platform users is available for test purpose only.
[LipidLynxX v0.4.12-beta preview release for Windows 10.](https://github.com/SysMedOs/LipidLynxX/releases/tag/v0.4.12-beta)

For macOS users, a installation pack is under development and will be ready in approximately end of May 2020.
If you really want to have an early access to the exe version, please contact us by email.

### Additional notice to developers

Since the code is still changing rapidly, the definitions of API and documentations in the source code may not be updated accordingly.
We kindly ask, if you have any plans to use LipidLynxX API contact us first, or follow this repository to get timely notifications when new changes are introduced.

### Key Features

- Optimized for manual interpretation and computer processing
- Suitable for both unmodified lipids and modified lipids
- Unified modification controlled vocabularies
- Unified position specific annotations
- Cross level match based on shared levels
- Extract key information from LipidLynxX ID
- Strictly controlled format using JSON schema
- Easy to use Graphic User Interface
- API access for professional users
- Command line tools for professional users


### Main Modules

- **LipidLynxX Converter**

  - Convert different abbreviations to uniformed LipidLynxX ID

- **LipidLynxX Equalizer**

  - Cross link different level of LipidLynxX ID on selected level

### LipidLynxX Nomenclature

- LipidLynxX levels

  - Lipid level:
    - **B**: Bulk
    - **D**: Discrete
    - **S**: sn Specific
  - Modification levels:

    - 0 : no modification
    - 1 : mass shift
    - 2 : element shift
    - 3 : number and type of modification
    - 4 : modification position information
    - 5 : additional information (e.g. R-/S-)

  - Double bond levels:

    - .0 : no information of double bond position (.0 should always be skipped, e.g. B0.0 -> B0)
    - .1 : double bond position information given
    - .2 : cis- / trans- information of all C=C bond

- LipidLynxX level matrix

  - The combinations of 3 sub-levels result in a matrix of LipidLynxX levels
    - e.g. B2 , D4, and S4.2

    | Mod  | DB   |      | Bulk   | Discrete |          |          | sn Specific |          |          |
    | ---- | ---- | ---- | ------ | -------- | -------- | -------- | ----------- | -------- | -------- |
    | 0    |      |      | **B**  | **D**    |          |          | **S**       |          |          |
    |      | .1   |      |        |          | **D0.1** |          |             | **S0.1** |          |
    |      |      | .2   |        |          |          | **D0.2** |             |          | **S0.2** |
    | 1.   |      |      | **B1** | **D1**   |          |          | **S1**      |          |          |
    |      | .1   |      |        |          | **D1.1** |          |             | **S1.1** |          |
    |      |      | .2   |        |          |          | **D1.2** |             |          | **S1.2** |
    | 2.   |      |      | **B2** | **D2**   |          |          | **S2**      |          |          |
    |      | .1   |      |        |          | **D2.1** |          |             | **S2.1** |          |
    |      |      | .2   |        |          |          | **D2.2** |             |          | **S2.2** |
    | 3.   |      |      | **B3** | **D3**   |          |          | **S3**      |          |          |
    |      | .1   |      |        |          | **D3.1** |          |             | **S3.1** |          |
    |      |      | .2   |        |          |          | **D3.2** |             |          | **S3.2** |
    | 4.   |      |      |        | **D4**   |          |          | **S4**      |          |          |
    |      | .1   |      |        |          | **D4.1** |          |             | **S4.1** |          |
    |      |      | .2   |        |          |          | **D4.2** |             |          | **S4.2** |
    | 5.   |      |      |        | **D5**   |          |          | **S5**      |          |          |
    |      | .1   |      |        |          | **D5.1** |          |             | **S5.1** |          |
    |      |      | .2   |        |          |          | **D5.2** |             |          | **S5.2** |

  - Example

    ![LipidLynx_01_Home](lynx/static/images/levels_mod_full.png)

- Currently supported modification controlled vocabularies
    ![LipidLynx_01_Home](doc/images/nomenclature_cv.png)

- Some examples of LipidLynx abbreviations:

  - Fatty acids

    - FA18:0
    - O-16:0
    - P-18:0
    - 20:4/<2OH,oxo>
    - 20:4/<{5Z,9E,11Z,14Z},OH{8S}>
    - 20:4/<{5Z,9E,12E,15E},2OH{8S,11R},oxo{14}>

  - Phospholipids
  
    - PC(O-16:0/18:1)
    - PE(P-16:0_18:1)
    - PC(16:0/20:4/<2OH,oxo>)
    - PE(16:0/20:4/<{5,9,12,15},2OH{8,11},oxo{14}>)

## Instructions

### Sample files:

- Test input file: `LipidLynxX/doc/sample_data/input`
- Test output file: `LipidLynxX/doc/sample_data/output`

### How to install and use LipidLynxX

Please find our user guide in folder `doc`.
-  [User Guide in PDF format](doc/LipidLynxX_UserGuide.pdf)
-  [User Guide in Markdown format](doc/LipidLynxX_UserGuide.md)

### Information for developers

- LipidLynxX is configured to use [travis-ci](https://travis-ci.com) and GitHub Actions with `py.test` to test
cross-platform compatibility on Linux, macOS and  Windows.

- Current status of the master branch 

    [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
    [![Codacy Badge](https://api.codacy.com/project/badge/Grade/f2180cda82034653ba57eed4473ed135)](https://app.codacy.com/gh/SysMedOs/LipidLynxX?utm_source=github.com&utm_medium=referral&utm_content=SysMedOs/LipidLynxX&utm_campaign=Badge_Grade_Dashboard) 
    ![GitHub commits since latest release](https://img.shields.io/github/commits-since/SysMedOs/LipidLynxX/v0.4.12-beta.svg)
    
    [![Travis-CI Build Status](https://travis-ci.com/SysMedOs/LipidLynxX.svg?branch=master)](https://travis-ci.com/SysMedOs/LipidLynxX)
    ![GitHub Actions Python application](https://github.com/SysMedOs/LipidLynxX/workflows/Python%20application/badge.svg)
    
- You can also use py.test to test LipidLynxX in your python environment, all test files can be found in `./test` folder.

### Errors/bugs

In case you experienced any problems with running LipidLynxX, 
please report an issue in the [issue tracker](https://github.com/SysMedOs/LipidLynxX/issues) or contact us.

### Screenshots

- **GUI**
    ![LipidLynx_02_Converter](doc/images/LipidLynxX_01_Converter_text_output.png)
- **API**

    Examples:
    ```bash
    curl http://127.0.0.1:5000/lynx/api/0.1/converter/str/ -d 'data="PLPC"' -X GET
    curl http://127.0.0.1:5000/lynx/api/0.1/converter/list/ -d 'data=["PAPE", "PE 36:4"]' -X GET
    curl http://127.0.0.1:5000/lynx/api/0.1/converter/dict/ -d 'data={"Sample1":["PAPC","PS 16:0/18:2(9Z,12Z)"],"SAMPLE2":["PC 16:0_20:4","PS 16:0_18:2"]}' -X GET
    curl http://127.0.0.1:5000/lynx/api/0.1/equalizer/ -d 'data={"Sample1":["PAPC","PS 16:0/18:2(9Z,12Z)", "DPPE"],"SAMPLE2":["PC 16:0_20:4","PS 16:0_18:2", "PG 18:0_18:3"]}&level=D0' -X GET
    ```
     ![LipidLynxX_API](doc/images/LipidLynxX_API.png)
    
    
- **Terminal Tools**
    
    - LipidLynxX Converter
    
    ```bash
    python LynxConverter.py -i doc/sample_data/input/LipidLynxX_test.xlsx -o doc/sample_data/output/LipidLynxX_test_converter_out.xlsx
    ```

    - LipidLynxX Equalizer

    ```bash
    python LynxEqualizer.py -l "B0,D0,D1" -i doc/sample_data/input/LipidLynxX_test.csv -o doc/sample_data/output/LipidLynxX_test_equalizer_out.xlsx
    ```

### License

- LipidLynxX is Dual-licensed

  - For academic and non-commercial use: `GPLv2 License`:

    - [The GNU General Public License version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

  - For commercial use: please contact the develop team by email.

- Please cite our publication in an appropriate form.

  - LipidLynxX preprint on `bioRxiv.org`

    - Zhixu Ni, Maria Fedorova.
        "LipidLynxX: lipid annotations converter for large scale lipidomics and epilipidomics datasets"

      - DOI: [10.1101/2020.04.09.033894](https://www.biorxiv.org/content/10.1101/2020.04.09.033894v1)

  - LipidLynx is based on the previous project [epiLION](https://github.com/SysMedOs/epiLION)

    - Ni, Zhixu, Laura Goracci, Gabriele Cruciani, and Maria Fedorova.
        "Computational solutions in redox lipidomics–Current strategies and future perspectives."
        Free Radical Biology and Medicine (2019).
      - DOI: [10.1016/j.freeradbiomed.2019.04.027](https://www.sciencedirect.com/science/article/pii/S0891584919303466)

### Report issues

- Report any issues here: <https://github.com/SysMedOs/LipidLynxX/issues>

### Fundings

We acknowledge all projects that supports the development of LipidLynxX:

- BMBF - Federal Ministry of Education and Research Germany:

    <https://www.bmbf.de/en/>

- e:Med Systems Medicine Network:

    <http://www.sys-med.de/en/>

- SysMedOS Project :

    <https://home.uni-leipzig.de/fedorova/sysmedos/>
