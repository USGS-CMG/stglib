Installation
************

Ensure you have a working Anaconda or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ installation. Both approaches below require this.

The easy way
============

Create a Pangeo CHS account.

Create an stglib environment:

``conda create -n stglib``

Activate the environment:

``conda activate stglib``

Install stglib:

``conda install -c conda-forge stglib``

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
