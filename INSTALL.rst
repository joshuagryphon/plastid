Plastid installation
====================

Bioconda
--------
.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
   :target: http://bioconda.github.io/recipes/plastid/README.html
   :alt: install with bioconda

``Bioconda`` is a channel for the conda package manager with a focus on
bioinformatics software. Once you have ``Bioconda`` installed, installing
``plastid`` is as easy as running

    $ conda create -n plastid plastid
    $ source activate plastid

This will install all of the necesary dependencies for ``plastid`` in an
isolated environment.


Requirements
------------


Non-python
..........

Compiling ``plastid`` requires a standard C build system (e.g. GCC or clang,
plus the standard C library, and linking tools) as well as ``zlib``. The good
news is that these are already required by ``plastid``'s dependencies, such as
``numpy`` and ``pysam``, so you probably already have them.


Python
......

Plastid requires Python 2.7 or 3.3 or greater, as well as a number of Python
packages. All but three of these (``numpy``, ``pysam``, and ``Cython``) --
will be automatically handled by Python's package management systems (``pip``
or ``easy_install`` if you are using them). ``numpy``, ``pysam`` and ``Cython``
are required for ``plastid``'s setup script to run.


Compiling & installing
----------------------

 #. First, make sure you have GCC or clang, zlib and its headers, and the
    Python headers. On an Ubuntu system, you could get these by typing::

        $ sudo apt-get install build-essential zlib1g zlib1g-dev python-dev

 #. Install ``numpy``, ``pysam``, and ``Cython``::

        $ sudo pip install numpy pysam cython

    or, for a single user install::

        $ pip install --user numpy pysam cython

 #. Compile & install ``plastid``::

        $ python setup.py build
        $ sudo python setup.py install

    or, for a single user install::

        $ python setup.py build
        $ python setup.py install --user



Building the documentation
--------------------------

Building the documentation requires plastid to be compiled so that its Cython
extensions are importable. In addition,  ``sphinx``, ``sphinxcontrib-bibtex``,
``sphinxcontrib-argdoc``, and ``numpydoc`` are required. Install these::

    $ sudo pip install sphinx sphinxcontrib-bibtex sphinxcontrib-argdoc numpydoc mock sphinx_rtd_theme

Then make the documentation::

    $ cd docs
    $ make html



Links
-----
  - `plastid <plastid.readthedocs.org>`_
  - `Bioconda <bioconda.github.io>`_
  - `GCC <gcc.gnu.org>`_
  - `clang <clang.llvm.org>`_
  - `zlib <www.zlib.net>`_
  - `Cython <cython.org>`_
  - `numpy <www.numpy.org>`_
  - `pysam <pysam.readthedocs.org>`_
  - `Sphinx <www.sphinx-doc.org>`_
  - `sphinxcontrib-bibtex <sphinxcontrib-bibtex.readthedocs.org>`_
  - `sphinxcontrib-argdoc <sphinxcontrib-argdoc.readthedocs.org>`_
  - `numpydoc <docs.scipy.org/doc/numpy-1.10.0/reference>`_
