# Pipeline for aligning Time series images using Metashape and arosics

Lorem Ipsum

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

## Introduction



## Features



## Installation

### Dowloading python, R, RStudio, Anaconda
- Python : https://www.python.org/downloads/
- Anaconda : https://www.anaconda.com/download/
- R and RStudio : https://posit.co/download/rstudio-desktop/

### Env python  

Open an anaconda prompt and type the following commands :

- conda env create -n env_name --file='path_to_cloned_repo/scripts_TS_arosics/environment.yaml'
- conda activate env_name
- pip install 'path_to_cloned_repo/scripts_TS_arosics/Metashape-2.1.0-cp37.cp38.cp39.cp310.cp311-none-win_amd64.whl'


### Env R 

Install reticulate package : type install.packages('reticulate') in the console

### How to execute code in R using the correct python environment
Two solutions :
- Copy the path of the environment you just created in the beginning of the script "exec_python_recalage.R"
(you can access it by using the command 'conda info --envs' in an anaconda prompt). The python env will be used only for this script.
- If you're using RStudio, go to Tools -> Global Options -> Python -> Select... -> Conda Environments, and select the env you just created. 
Warning : It will become the default python interpreter for all your R scripts 

### Activate Metashape key 

Add the line "scan.License().activate('your_license_key')" at the beginning of 2023-07-13_TimeSIFT_scripts_auto.py, right after the imports


## Usage



## License


