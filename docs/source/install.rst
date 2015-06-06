Installation
============

From PyPi
---------

Stable versions of :py:data:`yeti` can be fetched from `PyPi`_ using `Pip`_.

Simply type from the terminal::

    $ sudo pip install yeti


or, for a single-user install::

    $ pip install --user yeti


`pip`_ will fetch any uninstalled :ref:`dependencies <install_dependencies>`,
**unless** you are installing in a `Virtualenv`_. In this case you must
manually install dependencies first. For instructions, see in :ref:`install_fails_virtualenv`.


Command-line scripts will be installed wherever your system configuration dictates.
Typically the install path for command line scripts for users appears in
``~/bin`` or ``~/.local/bin``. For system-wide installs, they would be
in ``/usr/local/bin``. Make sure the appropriate location is in your ``PATH`` by
adding to your ``.bashrc``, or ``.profile``::

	export PATH=~/bin:~/.local.bin:/usr/local/bin:$PATH


Development versions
--------------------

Development versions can be fetched from `our github link`_::

    $ git clone git://`our github link`_



.. _install_dependencies :

Dependencies
------------
The following packages and their dependencies are required:

- Python
    - `Python`_     >= 2.7
    - `NumPy`_      >= 1.5
    - `SciPy`_      >= 0.12
    - `matplotlib`_ >= 1.3
    - `Pandas`_     >= 0.14.1
    - `Pysam`_      >= 0.7.7
    - `Biopython`_  >= 1.5
- Non-Python
    - `SAMtools`_   >= 0.1.19


The dependencies are optional, recommended, or required for specific functions:

- Python
	- `nose`_ for running the test suite in :py:mod:`yeti.test`
- Non-Python
	- `bowtie`_ (not Bowtie 2), for :py:mod:`~yeti.bin.crossmap`
	- `Jim Kent's utilities`_ for converting BED to BigBed files


