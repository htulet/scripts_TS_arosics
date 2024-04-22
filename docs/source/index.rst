.. TS_arosics_alignement documentation master file, created by
   sphinx-quickstart on Fri Apr 19 09:49:34 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TS_arosics_alignement's documentation!
=================================================
**TS_arosics_alignement** is a Python library that allows multi-epoch registration of images from the same area. It enables the creation and alignement of orthomosaics based on drone photos 
using the `Time-SIFT <https://arxiv.org/abs/1807.09700>`_ approach, and the use of `arosics <https://danschef.git-pages.gfz-potsdam.de/arosics/doc/>`_ to co-register the created orthomosaics 
and a reference image. Its tools can be called both in python and in R.
An `Agisoft Metashape <https://www.agisoft.com/>`_ license is necessary to call the Time-SIFT process.

Check out the :doc:`usage` section for further information, including how to
:ref:`install <installation>` the project.

.. note::

   This project is under active development.

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Contents
--------

.. toctree::

   usage