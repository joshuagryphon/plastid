Welcome to `plastid`!
====================

We're in beta!
--------------
Welcome! ``plastid`` is still in beta. This means that the package
is functional, but we're adding speed/memory optimizations, and considering
moving some things around inside the package so that they either make more
sense or take lest time to get to. We welcome your feedback
in this process, should you have any.

For documentation, see `our home page <http://plastid.readthedocs.org/en/latest/>`_
on `ReadtheDocs.org <http://readthedocs.org>`_.

To run the tests, download the :downoad:`test dataset <https://www.dropbox.com/s/h17go7tnas4hpby/plastid_test_data.tar.bz2?dl=0>` and unpack
it into ``plastid/test``.



Introduction
------------

``plastid`` is a Python library for genomic analysis -- in particular,
high-throughput sequencing data -- with an emphasis on simplicity for
users. It was written by Joshua Dunn in `Jonathan Weissman's lab <http://weissmanlab.ucsf.edu>`_
at `UCSF <http://ucsf.edu>`_,  initially for analysis of
ribosome profiling and RNA-seq data. Versions of it have been used
in several publications.

``plastid``'s intended audience includes computational and traditional biologists,
software developers, and even those who are new to sequencing analysis. It is
released under the :doc:`BSD 3-Clause license <license>`.

This package provides:

  #. A set of scripts that implement common sequencing
     analyses

  #. A set of classes that create a simple,
     intuitive interfaces to genomic features,
     read alignments, and quantitative data. These objects readily
     interace with existing scientific tools, like the SciPy stack.

  #. Script writing tools that make it easy to use the objects
     implemented in ``plastid``.

  #. Extensive documentation, both in source code and at readthedocs


Installation
------------
We're in development, so this takes a few extra steps:

    1. Make sure you have numpy and cython installed::

        pip install numpy cython

    2. Clone this repo::
        
        git clone https://github.com/joshuagryphon/plastid.git

    3. Install dependencies::

        cd plastid && pip install -r requirements.txt

    4. Build extensions::

        python setup.py build_ext --inplace

    5. Install beta::

        python setup.py develop --user
