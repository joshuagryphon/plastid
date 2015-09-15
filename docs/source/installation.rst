Installation
============


 .. commented for now
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


Get the Beta release
--------------------
We're in beta right now, so we haven't put anything on `PyPI`_ yet. To get
the development version, clone it from `our github repository <plastid_repo>`_.
From the terminal:

 .. code-block:: shell

    # get the source
    $ cd <some folder where you want to put plastid>
    $ git clone git://github.com/joshuagryphon/plastid.git

    # Install prereqs. Due to install peculiarities of numpy & scipy,
    # this takes a few more steps than we'd like:

    # get cython & numpy first
    $ pip install --user --upgrade cython numpy

    # then get scipy & the rest of the requirements
    $ pip install --user -r plastid/requirements.txt

    # build the cython extensions
    $ cd plastid
    $ python setup.py build_ext --inplace

    # install for a single user
    $ python setup.py develop --user



Non-Python Dependencies
-----------------------

The following are required or recommended for specific functions:

- Non-Python
   - `bowtie`_ (not Bowtie 2), for :py:mod:`~yeti.bin.crossmap`
   - `Jim Kent's utilities`_ for converting BED to BigBed files


