Installation
************

Obtain stglib by cloning the GitHub repo. Change to a directory where you'd like stglib to live and type:

``git clone https://github.com/dnowacki-usgs/stglib.git``

Ensure you have a working Anaconda or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ installation. After ``cd``ing to the directory containing stglib, type:

``conda env create -n stglib --file requirements-py37.yml``

This will create a Conda environment with the requirements for stglib installed. Activate the stglib environment by typing:

``conda activate stglib``

Then type:

``pip install -e . --no-deps``

This will install stglib and you can start processing data!
