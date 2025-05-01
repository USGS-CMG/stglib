Installation
************

The easy (local machine) way
============================

We recommend managing your Python packages using Miniforge.

Miniforge
----------

Install `Miniforge <https://conda-forge.org/download/>`_ by selecting the download appropriate for your platform.

Install and activate an stglib environment:

::

  conda create -n stglib stglib
  conda activate stglib

To update stglib

::

   conda activate stglib
   conda update stglib

Now you can start processing data!

If you want to contribute to stglib development
===============================================

Set up Miniforge as above. Obtain stglib by cloning the GitHub repo. Change to a directory where you'd like stglib to live and type:

``git clone https://code.usgs.gov/cmgp/stglib.git``

After ``cd``\ing to the stglib directory (``cd stglib``), type:

``conda env create -n stglib --file requirements.yml``

This will create a Conda environment with the requirements for stglib installed. Activate the stglib environment by typing:

``conda activate stglib``

Then type:

``pip install -e . --no-deps``

This will create an editable stglib installation so you can make changes to the codebase. Get the latest changes to stglib by running ``git pull``.

When new run scripts are added to stglib (for example, if support for a new instrument has been added), a ``git pull`` will not add them to your path. In this case you also need to re-run the ``pip install`` line above to install any new run scripts.

Missing CF standard names
=========================

If you get errors regarding undefined standard_names in your compliance-checker report, you may need to update the CF Standard Name Table. This is not done automatically, even when upgrading stglib. To update, run the following from the command line in an environment with stglib activated.

::

  conda update compliance-checker
  compliance-checker -d latest
