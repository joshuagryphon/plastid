Installation
============

From PyPi
---------
Stable versions of :py:data:`plastid` can be fetched from `PyPi`_ using `Pip`_.
Due to some quirks in Python packaging `numpy`_ and `pysam`_ must be installed
first:

Simply type from the terminal:

 .. code-block:: shell

    $ sudo pip install numpy pysam
    $ sudo pip install plastid


or, for a single-user install:

 .. code-block:: shell

    $ pip install --user numpy pysam
    $ pip install --user plastid


Command-line scripts will be installed wherever your system configuration dictates.
Typically the install path for command line scripts for users appears in
``~/bin`` or ``~/.local/bin``. For system-wide installs, they would be
in ``/usr/local/bin``. Make sure the appropriate location is in your ``PATH`` by
adding to your ``.bashrc``, or ``.profile``:

 .. code-block:: shell

   export PATH=~/bin:~/.local.bin:/usr/local/bin:$PATH


Development versions
--------------------
To fetch the latest development versions, clone it from
`our github repository <plastid_repo>`_. From the terminal:

 .. code-block:: shell

    # get the source
    $ cd <some folder where you want to put plastid>
    $ git clone git://github.com/joshuagryphon/plastid.git

    # Install prereqs. Due to install peculiarities of numpy & scipy,
    # this takes a few more steps than we'd like:

    # get cython & numpy first
    $ pip install --user --upgrade cython numpy

    # build the cython extensions
    $ python setup.py develop --user



Non-Python Dependencies
-----------------------

The following are required or recommended for specific functions:

   - `bowtie`_ (not Bowtie 2), for :py:mod:`~plastid.bin.crossmap`

   - `Jim Kent's utilities`_ for converting BED to BigBed files


