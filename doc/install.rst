Installation
************

The easy (CHS) way
==================

Create a Pangeo CHS account.

Log in to CHS Pangeo JupyterHub, and edit/create your ``~/.condarc`` file:

  channels:
    - conda-forge
    - defaults
  channel_priority: strict
  envs_dirs:
    - /home/jovyan/my-conda-envs

Then create the stglib environment, install stglib, and activate the environment:

  conda create -n stglib stglib ipykernel
  conda activate stglib

Now you can start processing data!

The easy (local machine) way
==================

Ensure you have a working Anaconda or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ installation.

Set your channel priority to conda-forge:

  conda config --add channels conda-forge
  conda config --set channel_priority strict

Install and activate an stglib environment:

  conda create -n stglib stglib
  conda activate stglib

Now you can start processing data!

If you want to contribute to stglib development
===============================================

Obtain stglib by cloning the GitHub repo. Change to a directory where you'd like stglib to live and type:

``git clone https://github.com/dnowacki-usgs/stglib.git``

After ``cd``ing to the directory containing stglib, type:

``conda env create -n stglib --file requirements-py37.yml``

This will create a Conda environment with the requirements for stglib installed. Activate the stglib environment by typing:

``conda activate stglib``

Then type:

``pip install -e . --no-deps``

This will create an editable stglib installation so you can make changes to the codebase.
