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
    $ git clone git://github.com/joshuagryphon/plastid.git

    # Do to a quirk in Python setup scripts, numpy,
    # and pysam must must be installed first:
    $ pip install --user --upgrade numpy pysam

    # Install in develop mode
    $ cd plastid
    $ python setup.py develop --user


Non-Python Dependencies
-----------------------

The following are required or recommended for specific functions:

   - `bowtie`_ (not Bowtie 2), for :py:mod:`~plastid.bin.crossmap`

   - `Jim Kent's utilities`_ for converting BED to BigBed files


