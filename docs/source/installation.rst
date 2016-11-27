Installation
============

.. contents::
   :local:
 


From PyPi (recommended)
-----------------------

Install package
...............

Stable versions of :data:`plastid` can be fetched from `PyPi`_ using `Pip`_.
Due to some quirks in Python packaging `Cython`_, `numpy`_ , `pysam`_ must be
installed first:

Simply type from the terminal:

.. code-block:: shell

   $ sudo pip install numpy pysam cython
   $ sudo pip install plastid

or, for a single-user install:

.. code-block:: shell

   $ pip install --user numpy pysam cython
   $ pip install --user plastid

Test your installation within Python:

.. code-block:: python

   >>> import plastid

If installation fails with a message that an object in `numpy`_ or `pysam`_ is
the wrong size, the included C files may need to be regenerated.

To do so, first make sure all the dependencies are installed. Then type:

.. code-block:: shell

   $ pip install --verbose --user --install-option="--recythonize" plastid

And then re-test the installation. If installation continues to fail, please see
:ref:`faq-install-fails` for common errors or `our issue tracker`_ to report a
new one.


Set ``PATH`` variable
.....................
Command-line scripts will be installed wherever system configuration dictates. On OSX and many varities of linux, the install path for a single-user install is ``~/bin`` or ``~/.local/bin``. For system-wide installs, the path is typically ``/usr/local/bin``. Make sure the appropriate location is in your ``PATH`` by adding to your ``.bashrc``, ``.bash_profile``, or ``.profile``:

.. code-block:: shell

    export PATH=~/bin:~/.local.bin:/usr/local/bin:$PATH

Also, type the line above in any open terminal (or login and out again) to apply the changes.


.. _install-inside-venv:

Inside a `virtualenv`
---------------------

Often users or systems administrators need to install multiple versions of the
same package for different scientific purposes. To do so they use *sandboxes*
that insulate packages from each other.

The easiest wsay to install :data:`Plastid` inside a sandbox is to use
`virtualenv`_:

.. code-block:: shell

   # install virtualenv if you don't have it.
   # use either "sudo" or "--user", not both.

   # Use this line for a system-wide install
   $ sudo pip install virtualenv

   # or, use this line for single user install
   $ pip install --user virtualenv

   # With virtualenv installed, create & activate vanilla environment
   # when prompted, do NOT give the virtualenv access to system packages

   # create
   $ virtualenv ~/some/path/to/venv

   # activate
   $ source ~/some/path/to/venv/bin/activate

   # Fresh install of plastid.
   # Note- no use of `sudo` here. It confuses the virtualenv
   (venv) $ pip install numpy pysam cython
   (venv) $ pip install plastid

   # test
   (venv) $ python -c "from plastid import *"



Development versions
--------------------
To fetch the latest development versions, clone it from `our github repository <plastid_repo>`_. From the terminal:

.. code-block:: shell

   # get the source
   $ git clone git://github.com/joshuagryphon/plastid.git

   # Do to a quirk in Python setup scripts, numpy,
   # and pysam must must be installed first:
   $ pip install --user --upgrade numpy pysam

   # Install in develop mode
   # Use `--recythonize` flag to link code against your
   # versions of numpy and pysam, if they are different
   # from ours
   $ cd plastid
   $ python setup.py develop --user --recythonize


Non-Python Dependencies
-----------------------

Plastid has a number of non-Python dependencies:

 - A full build system for C compiling (e.g. `GCC <gcc.gnu.org>`_ or `clang <clang.llvm.org>`_) 
 - `zlib <www.zlib.net>`_, including its headers


The following are not required for full functionality, but are recommended for specific functions:

 - `bowtie`_ (not `bowtie 2`_) for use in  :py:mod:`~plastid.bin.crossmap`
 - `Jim Kent's utilities`_ for converting BED to BigBed files
 - `The FASTX toolkit <http://hannonlab.cshl.edu/fastx_toolkit/>`_   



Troubleshooting
---------------

:data:`plastid` installs fairly easily in most Linux and Macintosh setups.

One exception is under `Anaconda`_, which, depending upon your setup, can require 
custom work to get going. If you need to run :data:`plastid` inside a sandbox,
we strongly recommend using `virtualenv`_ rather than `Anaconda`_ for this
purpose. To do so, see :ref:`install_inside_venv`, above.

For other troubleshooting, please see our FAQ section on :ref:`installation <faq-run>`.
