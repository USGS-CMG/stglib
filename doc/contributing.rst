.. _contributing:

**********************
Contributing to stglib
**********************

.. contents:: Table of contents:
   :local:

.. note::

  We used `Contributing to xarray <http://xarray.pydata.org/en/stable/contributing.html>`_ as a guide,
  which in turn came from the `Pandas Contributing
  Guide <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_.

Where to start?
===============

We need your help.  All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

If you are brand new to open-source development, we recommend going
through the `stglib's GitHub "issues" tab <https://github.com/dnowacki-usgs/stglib/issues>`_
to find issues that interest you, discuss issues, to report new issues or propose new features.


.. _contributing.bug_reports:

Bug reports and enhancement requests
====================================

To indicate desire to work on an issue, post your intent and ideas in the issues tab.
Follow the guidance here to set up your environment in python.

We do not yet have a public mailing list.  USGS folks can ask questions and follow discussions in the
GS-CMHRP-CCH-TimeSeriesData team, however we encourage all to publicly post and discuss issues here
on github in `stglib's GitHub "issues" tab <https://github.com/dnowacki-usgs/stglib/issues>`_,
so that we can advance this package as a user community.

Bug reports are an important part of improving our package. Having a complete bug
report will allow others to reproduce the bug and provide insight into fixing. See
`this stackoverflow article <https://stackoverflow.com/help/mcve>`_ for tips on
writing a good bug report.

Trying the bug-producing code out on the *master* branch is often a worthwhile exercise
to confirm the bug still exists. It is also worth searching existing bug reports and
pull requests to see if the issue has already been reported and/or fixed.

Bug reports must:

#. Include a short, self-contained Python snippet reproducing the problem.
   You can format the code nicely by using `GitHub Flavored Markdown
   <http://github.github.com/github-flavored-markdown/>`_
#. Explain why the current behavior is wrong/not desired and what you expect instead.

The issue will then show up to the *stglib* community and be open to comments/ideas
from others.

.. _contributing.github:

Working with the code
=====================

Now that you have an issue you want to fix, enhancement to add, or documentation
to improve, you need to learn how to work with GitHub and the *stglib* code base.

.. _contributing.version_control:

Version control, Git, and GitHub
--------------------------------

To contribute, you will need to know git.
To the new user, working with Git is one of the more daunting aspects of contributing
to *stglib*.  It can very quickly become overwhelming, but sticking to the guidelines
below will help keep the process straightforward and mostly trouble free.  As always,
if you are having difficulties please feel free to ask for help.

The code is hosted on `GitHub <https://github.com/dnowacki-usgs/stglib>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <http://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning Git:

* the `GitHub help pages <http://help.github.com/>`_.
* the `NumPy's documentation <http://docs.scipy.org/doc/numpy/dev/index.html>`_.
  "Every single developer working on the project has their code reviewed, and we\'ve
  come to see it as friendly conversation from which we all learn and the overall code
  quality benefits. Therefore, please don’t let the review discourage you from contributing:
  its only aim is to improve the quality of project, not to criticize (we are, after all,
  very grateful for the time you’re donating!)."
* Matthew Brett's `Curious Coder\'s Guide To Git <https://matthew-brett.github.io/curious-git/>`_.

Getting started with Git
------------------------

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _contributing.forking:

Forking
-------

You will need your own fork to work on the code. Go to the `stglib project
page <https://github.com/dnowacki-usgs/stglib>>`_ and hit the ``Fork`` button. You will
want to clone your fork to your machine::

    git clone https://github.com/your-github-user-name/stglib.git
    cd stglib
    git remote add upstream https://github.com/dnowacki-usgs/stglib.git

This creates the directory `stglib` and connects your repository to
the upstream (main project) *stglib* repository.

.. _contributing.dev_env:

Creating a development environment
----------------------------------

To test out code changes, you'll need to build *stglib* from source, which
requires a Python environment. If you're making documentation changes, you can
skip to `contributing.documentation` but you won't be able to build the
documentation locally before pushing your changes.

.. _contributing.dev_python:

Creating a Python Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before starting any development, you'll need to create an isolated stglib
development environment:

- We recommend installing the IOOS conda environment per these instructions\:
  `Installing the IOOS Environment <http://ioos.github.io/notebooks_demos/other_resources/>`_
- Make sure your conda is up to date with the command (``conda update conda``)
- Make sure that you have cloned the repository
- ``cd`` to the *stglib* source directory (your fork, locally, on your own machine)
- install *stglib* per `Installation <https://stglib.readthedocs.io/en/latest/install.html>`

At this point you should be able to import *stglib* from your locally built version in
a python interpreter or in jupyter-notebook::

   $ python  # start an interpreter
   >>> import stglib
   >>> stglib.__version__
   '0.1.0+12.gd81f135'

The above procedure created a new environment, and did not touch any of your existing environments,
nor any existing Python installation.

To view your environments::

      conda info -e

To return to your root (or base) environment::

      conda deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

Creating a git branch
---------------------

You want your master branch to reflect only production-ready code, so create a
feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to *stglib*. You can have many "shiny-new-features" as individual branches
and switch in between those branches using the ``git checkout the-feature-branch`` command.

To update your shiny-new-feature branch, you need to retrieve the changes from the master branch::

    git fetch upstream
    git rebase upstream/master

Keep in mind, `upstream` refers to the original version of *stglib* at
`<https://github.com/dnowacki-usgs/stglib>`,
not to be confused with the term `origin`, which is your fork of *stglib* at
`<https://github.com/your-github-user-name/stglib.git>`.
The fetch and rebase commands will replay your commits (changes) on top of the latest *stglib* git master.
If this leads to merge conflicts, you must resolve these before submitting your pull
request.  If you have uncommitted changes that you are not ready to commit yet,
you will need to ``git stash`` them prior to updating.  ``git stash`` will effectively store your changes
and they can be reapplied with ``git stash pop`` after updating.

.. _contributing.documentation:

Contributing to the documentation
=================================

If you're not the developer type, contributing to the documentation is still of
huge value. You don't even have to be an expert on *stglib* to do so! In fact,
there are sections of the docs that are worse off after being written by
experts. If something in the docs doesn't make sense to you, updating the
relevant section after you figure it out is a great way to ensure it will help
the next person.

.. contents:: Documentation:
   :local:


About the *stglib* documentation
--------------------------------

The documentation is written in `reStructuredText <https://en.wikipedia.org/wiki/ReStructuredText>`__,
which is almost like writing in plain English, and built using `Sphinx <http://sphinx-doc.org/>`__. The
Sphinx Documentation has an excellent `introduction to reST
<http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__.
Review the Sphinx docs to perform more complex changes to the documentation as well.

Some other important things to know about the docs:

- The *stglib* documentation consists of two parts: the docstrings in the code
  itself and the docs in this folder ``stglib/doc/``.

  The docstrings are meant to provide a clear explanation of the usage of the
  individual functions, while the documentation in this folder consists of
  tutorial-like overviews per topic together with some other information
  (what's new, installation, etc).

- The docstrings follow the **Numpy Docstring Standard**, which is used widely
  in the Scientific Python community. This standard specifies the format of
  the different sections of the docstring. See `this document
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  for a detailed explanation, or look at some of the existing functions to
  extend it in a similar manner.

- `stglib` documentation is organized by instrument type

- There is an index for all the documentation called ``index.rst`` and if you make a new
  documentation file for some new instrument, be sure to include it in a ``toctree`` in ``index.rst``


How to build the *stglib* documentation
---------------------------------------

Requirements
~~~~~~~~~~~~
Follow the instructions on creating a development environment above, and to build the docs
you need to create a new environment with the environment file ``doc/doc-requirements.yml``.

.. code-block::

    # Create and activate the docs environment
    conda env create -f doc/doc-requirements.yml
    conda activate stglib-docs

    # Build and install stglib
    pip install -e .

Building the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to your local ``stglib/doc/`` directory in the console and run::

    make html

Then you can find the HTML output in the folder ``stglib/doc/build/html/``.

The first time you build the docs, it will take longer because it has to run
all the code examples and build all the generated docstring pages. In subsequent
evocations, sphinx will try to only build the pages that have been modified.

If you want to do a full clean build, do::

    make clean
    make html

.. _contributing.code:

Contributing to the code base
=============================

.. contents:: Code Base:
   :local:

Code standards
--------------

Writing good code is not just about what you write. It is also about *how* you
write it.

We expect any new code to be well documented, both in the code itself and for there
to be explanations and tutorials similar to what already exists in the ``doc/`` directory.

We expect new code to follow the structure of the existing code.

In addition, because a lot of people use our library, it is important that we
do not make sudden changes to the code that could have the potential to break
a lot of user code as a result, that is, we need it to be as *backwards compatible*
as possible to avoid mass breakages.

Code Formatting
~~~~~~~~~~~~~~~

stglib follows PEP8 conventions.

There are tools you can use to ensure a consistent code format:

- `Black <https://black.readthedocs.io/en/stable/>`_ for standardized code formatting
- `Flake8 <http://flake8.pycqa.org/en/latest/>`_ for general code quality
- `isort <https://github.com/timothycrosley/isort>`_ for standardized order in imports.
  See also `flake8-isort <https://github.com/gforcada/flake8-isort>`_.
- `mypy <http://mypy-lang.org/>`_ for static type checking on `type hints
  <https://docs.python.org/3/library/typing.html>`_

Integrated development environments also help with code formatting:

 - `spyder <https://www.spyder-ide.org/>`_  installed by ``conda install spyder``
 - `Atom <https://atom.io/>`_
 - `pycharm free community edition <https://www.jetbrains.com/pycharm/>`_ is very full featured and
   plays well with conda, however can be hard to learn
 - `vscode <https://code.visualstudio.com/>`_ is simpler, and may have issues with conda


Backwards Compatibility
~~~~~~~~~~~~~~~~~~~~~~~

Please try to maintain backward compatibility.  If you think breakage is
required, clearly state why as part of the pull request.  Also, be careful when changing
method signatures and add deprecation warnings where needed.

.. _contributing.ci:

Testing With Continuous Integration
-----------------------------------

We use continuous integration testing, which evaluates the code each time
code is ``pushed`` to github.

The *stglib* test suite consists of the files in ``stglib/tests/``, and are run automatically at
`Travis CI <https://travis-ci.org/dnowacki-usgs/stglib>`__,
a continuous integration service, once your pull request is submitted.

You may wish to run tests on your local branch before pushing to github or submitting the pull request.

There are several types of testing:

 - The simplest, and built into python, is
   `unittest <https://docs.python.org/2/library/unittest.html>`_.  ``test_stglib.py`` uses unittest.
   Try running the tests with the command ``python -m unittest discover stglib\tests``
 - `pytest <https://docs.pytest.org/en/latest/>`_ can be used for more complicated testing.
   ``test_puv_quick.py`` uses pytest.  You will need to install pytest (``conda install pytest``)
   before you can use it.  ``pytest`` can be run from within the ``tests`` directory.  pytest will run
   test written for pytest and for unittest.

A pull-request will be considered for merging when you have an all 'green' build. If any
tests are failing, then you will get a red 'X', where you can click through to see the
individual failed tests.

.. _contributing.tdd:


Test-driven development/code writing
------------------------------------

All tests should go into the ``tests`` subdirectory of the specific package.
This folder contains many current examples of tests, and we suggest looking to these for
inspiration.

`test-driven development (TDD) <http://en.wikipedia.org/wiki/Test-driven_development>`_:
This development process "relies on the repetition of a very short development cycle:
first the developer writes an (initially failing) automated test case that defines a desired
improvement or new function, then produces the minimum amount of code to pass that test."
So, before actually writing any code, you should write your tests.  Often the test can be
taken from the original GitHub issue.  However, it is always worth considering additional
use cases and writing corresponding tests.

*stglib* maintainers will ask that your code include tests when receiving a pull request.  Therefore,
it is worth getting in the habit of writing tests ahead of time so this is never an issue.

For more information about how to write tests, the xarray maintainers have `writing tests for xarray
<http://xarray.pydata.org/en/stable/contributing.html#test-driven-development-code-writing>`_

We will include more information here as stglib grows.

Contributing your changes to *stglib* (how to use git)
======================================================

Committing your code
--------------------

Keeping style fixes to a separate commit will make your pull request more readable.

Once you've made changes, you can see them by typing::

    git status

If you have created a new file, it is not being tracked by git. Add it by typing::

    git add path/to/file-to-be-added.py

Doing 'git status' again should give something like::

    # On branch shiny-new-feature
    #
    #       modified:   /relative/path/to/file-you-added.py
    #

The following defines how a commit message should be structured:

    * A subject line with `< 72` chars.
    * One blank line.
    * Optionally, a commit message body.

Please reference the relevant GitHub issues in your commit message using ``GH1234`` or
``#1234``.  Either style is fine, but the former is generally preferred.

Now you can commit your changes in your local repository::

    git commit -m

Squashing your commits
----------------------

*stglib* maintainers prefer that commits be ``squashed`` before a pull request is initiated.
This is difficult to do once commits are pushed to github.  This is explained `here
<https://stackoverflow.com/questions/5189560/squash-my-last-x-commits-together-using-git>`_.


Pushing your changes
--------------------

When you want your changes to appear publicly on your GitHub page, push your
forked feature branch's commits::

    git push origin shiny-new-feature

Here ``origin`` is the default name given to your remote repository on GitHub (your fork of stglib).
You can see the remote repositories::

    git remote -v

If you added the upstream repository as described above you will see something
like::

    origin  git@github.com:yourname/stglib.git (fetch)
    origin  git@github.com:yourname/stglib.git (push)
    upstream        git://github.com/pydata/stglib.git (fetch)
    upstream        git://github.com/pydata/stglib.git (push)

Now your code is on GitHub, but it is not yet a part of the *stglib* project.  For that to
happen, a pull request needs to be submitted on GitHub.

Review your code
----------------

When you're ready to ask for a code review, file a pull request. Before you do, once
again make sure that you have followed all the guidelines outlined in this document
regarding code style, tests, performance tests, and documentation. You should also
double check your branch changes against the branch it was based on:

#. Navigate to your repository on GitHub -- https://github.com/your-user-name/stglib
#. Click on ``Branches``
#. Click on the ``Compare`` button for your feature branch
#. Select the ``base`` and ``compare`` branches, if necessary. This will be ``master`` and
   ``shiny-new-feature``, respectively.

Finally, make the pull request
------------------------------

If everything looks good, you are ready to make a pull request.  A pull request is how
code from a local repository becomes available to the GitHub community and can be looked
at and eventually merged into the master version.  This pull request and its associated
changes will eventually be committed to the master branch and available in the next
release.  To submit a pull request:

#. Navigate to your repository on GitHub
#. Click on the ``Pull Request`` button
#. You can then click on ``Commits`` and ``Files Changed`` to make sure everything looks
   okay one last time
#. Write a description of your changes in the ``Preview Discussion`` tab
#. Click ``Send Pull Request``.

This request then goes to the repository maintainers, and they will review
the code. If you need to make more changes, you can make them in
your branch, add them to a new commit, push them to GitHub, and the pull request
will be automatically updated.  Pushing them to GitHub again is done by::

    git push origin shiny-new-feature

This will automatically update your pull request with the latest code and restart the
Travis Continuous Integration tests.


Delete your merged branch (optional)
------------------------------------

Once your feature branch is accepted into upstream, you'll probably want to get rid of
the branch. First, merge upstream master into your branch so git knows it is safe to
delete your branch::

    git fetch upstream
    git checkout master
    git merge upstream/master

Then you can do::

    git branch -d shiny-new-feature

Make sure you use a lower-case ``-d``, or else git won't warn you if your feature
branch has not actually been merged.

The branch will still exist on GitHub, so to delete it there do::

    git push origin --delete shiny-new-feature


PR checklist
------------

- **Properly comment and document your code.**
- **Test that the documentation builds correctly** by typing ``make html`` in the ``doc`` directory.
  This is not strictly necessary, but this may be easier than waiting for CI to catch a mistake.
- **Test your code**.

    - Write new tests if needed.
    - Test the code using or unittest.

- **Properly format your code**
- **Squash your commits**
- **Push your code and** `create a PR on GitHub <https://help.github.com/en/articles/creating-a-pull-request>`_.
- **Use a helpful title for your pull request** by summarizing the main contributions rather
  than using the latest commit message. If this addresses an `issue <https://github.com/dnowacki-usgs/stglib/issues>`_,
  please `reference it <https://help.github.com/en/articles/autolinked-references-and-urls>`_.
