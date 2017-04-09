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

   $ pip install --verbose --user --no-cache-dir --install-option="--recythonize" plastid

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

Inside a virtualenv
-------------------

Often users or systems administrators need to install multiple versions of the
same package for different scientific purposes. To do so they use *sandboxes*
that insulate packages from each other.

The easiest way to install :data:`Plastid` inside a sandbox is to use
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
   (venv) $ pip install --no-cache-dir plastid

   # test
   (venv) $ python -c "from plastid import *"



Development versions
--------------------
To fetch the latest development versions, clone it from `our github repository <plastid_repo>`_. From the terminal:

.. code-block:: shell

   # get the source
   $ git clone git://github.com/joshuagryphon/plastid.git

   # Do to a quirk in Python setup scripts (fixed in new versions of pip),
   # numpy, cython, and pysam must must be installed first:
   $ pip install --user --upgrade numpy pysam

   # Install in develop mode. Use `--recythonize` flag to regenerate
   # C files if necessary (e.g. after upgrading pysam)
   $ cd plastid
   $ pip install --install-option='--recythonize' --user -e .


Non-Python Dependencies
-----------------------

Plastid has a number of non-Python dependencies:

 - A full build system for C compiling (e.g. `GCC <gcc.gnu.org>`_ or `clang <clang.llvm.org>`_) 
 - `zlib <www.zlib.net>`_, including its headers


The following are not required for full functionality, but are required for specific functions or just plain useful:

 - `bowtie`_ (not `bowtie 2`_) for use in  :py:mod:`~plastid.bin.crossmap`
 - `Jim Kent's utilities`_ for converting BED to BigBed files
 - `The FASTX toolkit <http://hannonlab.cshl.edu/fastx_toolkit/>`_   



Troubleshooting
---------------

:data:`plastid` installs fairly easily in most Linux and Macintosh setups. If
you run into issues running or installing, please see our
FAQ section on :ref:`installation <faq-run>`
and then `our issue tracker`_ to see if anybody else has encountered your issue,
and if instructions already exist.

Frequently, problems can be solved by installing :data:`plastid` in a clean
environment. For instructions, see :ref:`install-inside-venv`, above.


Notes on `conda`_
.................

Numerous users have reported unexpected behaviors when installing
:data:`plastid` in `conda`_/`Anaconda`_ environments. These can appear
`at runtime <http://>`_,
or
`after installation <http://>`_.
Notably, :data:`plastid` is not the only package to experience this
incompatibility.

These issues can take a long time to resolve and frequently the solution is 
installation inside a clean environment that does not use `conda`_. We're
looking into this issue to improve our compatibility, but at present, `conda`_
is not supported. If you need to run :data:`plastid` inside a sandbox,
we strongly recommend using `virtualenv`_. To do so, see :ref:`install-inside-venv`, above.

For other troubleshooting, please see our FAQ section on :ref:`installation <faq-run>`.
