Installation
============

From PyPi
---------
Stable versions of :py:data:`yeti` can be fetched from `PyPi`_ using `Pip`_.

Simply type from the terminal:

 .. code-block:: shell

    $ sudo pip install yeti


or, for a single-user install:

 .. code-block:: shell

    $ pip install --user yeti


Command-line scripts will be installed wherever your system configuration dictates.
Typically the install path for command line scripts for users appears in
``~/bin`` or ``~/.local/bin``. For system-wide installs, they would be
in ``/usr/local/bin``. Make sure the appropriate location is in your ``PATH`` by
adding to your ``.bashrc``, or ``.profile``:

 .. code-block:: shell

   export PATH=~/bin:~/.local.bin:/usr/local/bin:$PATH


From PyPi in a `virtualenv`_
----------------------------
Due to `peculiarities in the setup requirements <https://github.com/numpy/numpy/issues/2434>`_
for `numpy`_ and `scipy`_, you must manually install these in your `virtualenv`_
**before** installing :data:`yeti`. See :ref:`install_fails_virtualenv` for instructions.


Development versions
--------------------
Development versions can be fetched from `our github link`_:

 .. code-block:: shell

    $ git clone git://



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


The following are required or recommended for specific functions:

- Python
   - `nose`_ for running the test suite in :py:mod:`yeti.test`
- Non-Python
   - `bowtie`_ (not Bowtie 2), for :py:mod:`~yeti.bin.crossmap`
   - `Jim Kent's utilities`_ for converting BED to BigBed files


