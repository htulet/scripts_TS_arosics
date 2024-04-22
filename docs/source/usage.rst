Usage
=====

.. _installation:

Installation
------------
Dowloading python, R, RStudio, Anaconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To use TS_arosics_alignement, first install Python (3.7 or superior), Anaconda and R :
- Python : `https://www.python.org/downloads/ <https://www.python.org/downloads/>`_
- Anaconda : `https://www.anaconda.com/download/ <https://www.anaconda.com/download/>`_ 
- R and RStudio : `https://posit.co/download/rstudio-desktop/ <https://posit.co/download/rstudio-desktop/>`_. *R and Rstudio are only useful if you want to execute the code using a R IDE.*

Python environment  
~~~~~~~~~~~~~~~~~~
You will then need to create a python environment with all the correct dependencies. Open an anaconda prompt and type the following commands

.. code-block:: console

   $ conda env create -n env_name --file='path_to_cloned_repo/scripts_TS_arosics/environment.yaml'
   $ conda activate env_name

To install the Metashape python API, download the wheel file corresponding to your operating system here : `https://www.agisoft.com/downloads/installer/ <https://www.agisoft.com/downloads/installer/>`_ 
Then, run the following command

.. code-block:: console

   (env_name) $ pip install 'path_to_downloaded_whl_file'

*Don't forget each time to replace the paths and the name of your environment by the actual path and names you want to use*

R Environment and how to link it to python (Optionnal)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you want to execute the code from an R IDE, you will need to link the R environment you use to the python environment you just created. 
Fist, install the reticulate package in your R environment using the Rstudio console

.. code-block:: console

   $ install.packages('reticulate')

Then, you can either:
- Copy the path of the environment you just created in the beginning of the script "exec_python_recalage.R"
(you can access it by using the command **conda info --envs** in an anaconda prompt). The python env will be used only for this script.
- In RStudio, you can go to Tools -> Global Options -> Python -> Select... -> Conda Environments, and select the env you just created. 
Warning : It will become the default python interpreter for all your R scripts 


Activation of the Metashape key 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Write down your license key at the beginning of the script 'TimeSIFT_scripts_auto.py', right after the imports, in the line "scan.License().activate('your_license_key')", then execute the script.

Generate orthomosaics from drone images using Time-SIFT and Metashape
---------------------------------------------------------------------

You can use the ``Time_SIFT_process()`` function

.. autofunction:: TimeSIFT_scripts_auto.Time_SIFT_process

Align two or more images together using arosics
-----------------------------------------------

You can use the ``arosics_chain.complete_arosics_process()`` function

.. autofunction:: arosics_chain.complete_arosics_process

