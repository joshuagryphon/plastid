Welcome to `plastid`!
=====================

For documentation, see `our home page
<http://plastid.readthedocs.io/en/latest/>`_ on `ReadtheDocs
<http://readthedocs.io>`_.

To run the tests, download the `test dataset
<https://www.dropbox.com/s/np3wlfvp6gx8tb8/2022-05-04.plastid-test-data.tar.bz2?dl=0>`_
and unpack it into ``plastid/test``.


Introduction
------------

``plastid`` is a Python library for genomic analysis -- in particular,
high-throughput sequencing data -- with an emphasis on simplicity for users. It
was written by Joshua Dunn in `Jonathan Weissman's lab
<http://weissmanlab.ucsf.edu>`_ at `UCSF <http://ucsf.edu>`_,  initially for
analysis of ribosome profiling and RNA-seq data. Versions of it have been used
in several publications.

``plastid`` intended audience includes computational and traditional
biologists, software developers, and even those who are new to sequencing
analysis. It is released under the BSD 3-Clause license.

This package provides:

#. A set of scripts that implement common sequencing analyses

#. A set of classes for exploratory data analysis. These provide simple
   and consistent interfaces for manipulating genomic features,
   read alignments, and quantitative data; and readily interface with
   existing scientific tools, like the SciPy stack.

#. Script writing tools that make it easy to use the objects implemented in
   ``plastid``.

#. Extensive documentation, both in source code and at `our home page
   <http://plastid.readthedocs.io/en/latest/>`_ on `ReadtheDocs
   <http://readthedocs.io>`_.


Installation
------------

Bioconda
........

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
   :target: http://bioconda.github.io/recipes/plastid/README.html
   :alt: install with bioconda

``Bioconda`` is a channel for the conda package manager with a focus on
bioinformatics software. Once you have ``Bioconda`` installed, installing
``plastid`` is as easy as running::

    $ conda create -n plastid plastid
    $ source activate plastid

This will install all of the necessary dependencies for ``plastid`` in an
isolated environment.

PyPI
....

``plastid`` can be installed directly from PyPI

    $ pip install numpy pysam cython
    $ pip install plastid

If you get any runtime warnings about numpy versions having changed, or about
a missing module in Pysam, or about some object being the wrong size, try
regenerating the included C source files from the original Cython code. To
do this type::

    $ pip install --upgrade  --install-option='--recythonize' plastid


Running the tests
-----------------

- **NOTE**: to run the entire test suite you'll first need to download our `test
  dataset`_ and unpack it into `plastid/test/data`.

We use ``nose`` as our test runner, and test under different versions of Python
using ``tox``. To completely control the environment (e.g. compilers et c), we
recommend running the tests inside the Docker container, which contains 
large data files needed for the tests that aren't packaged with ``plastid`` by
default:

.. code-block:: shell

   # build & run the Docker image from within the project folder
   $ docker build --pull -t plastid .
   $ docker run -it --rm plastid

   # inside the container, run the tests over all default configurations
   root@plastid $ tox


Our ``tox`` config lets developers run subsets of tests rather than the full
suite.  All positional arguments are passed through to ``nosetests``

.. code-block:: shell

   # run all tests within the plastid.test.unit subpackage
   root@plastid $ tox plastid.test.unit

   # run tests in two files
   root@plastid $ tox plastid.test.unit.genomics.readers.test_bed plastid.test.unit.util.io.test_binary

By default, ``tox`` recompiles all C extensions before running the tests. This
can be slow. To avoid doing that, set the environment variable
`PLASTID_NOREBUILD` to `true`:

.. code-block:: shell

   # run unit tests without rebuilding the C extensions
   root@plastid $ env PLASTID_NOREBUILD=true tox plastid.test.unit

Finally, if you only want to test in some, not all environments, you can do so
with typical ``tox`` syntax:

.. code-block:: shell

   # list available test environments
   root@plastid $ tox -l
   py36-pinned
   py36-latest
   py39-latest

   # run only in 2 selected environments
   root@plastid $ tox -e py36-pinned,py39-latest plastid.test.unit



Links & help
------------

- `Documentation <http://plastid.readthedocs.io>`_

- `Our github repo <https://github.com/joshuagryphon/plastid>`_

