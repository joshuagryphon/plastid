Welcome to `plastid`!
=====================

For documentation, see `our home page
<http://plastid.readthedocs.io/en/latest/>`_ on `ReadtheDocs
<http://readthedocs.io>`_.

To run the tests, download the `test dataset
<https://www.dropbox.com/s/h17go7tnas4hpby/plastid_test_data.tar.bz2?dl=0>`_
and unpack it into ``plastid/test``.


Introduction
------------

``plastid`` is a Python library for genomic analysis -- in particular,
high-throughput sequencing data -- with an emphasis on simplicity for users. It
was written by Joshua Dunn in `Jonathan Weissman's lab
<http://weissmanlab.ucsf.edu>`_ at `UCSF <http://ucsf.edu>`_,  initially for
analysis of ribosome profiling and RNA-seq data. Versions of it have been used
in several publications.

``plastid``'s intended audience includes computational and traditional
biologists, software developers, and even those who are new to sequencing
analysis. It is released under the BSD 3-Clause license.

This package provides:

  #. A set of scripts that implement common sequencing analyses

  #. A set of classes for exploratory data anlysis. These provide simple
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

``plastid`` can be installed directly from PyPI, but requires numpy, pysam,
and cython to be installed first i.e.::

    $ pip install numpy pysam cython
    $ pip install plastid

If you get any runtime warnings about numpy versions having changed, or about
a missing module in Pysam, or about some object being the wrong size, try
regenerating the included C source files from the original Cython code. To
do this type::

    $ pip install --upgrade  --install-option='--recythonize' plastid


Links & help
------------

  - `Documentation <http://plastid.readthedocs.io>`_

  - `Our github repo <https://github.com/joshuagryphon/plastid>`_

  - Subscribe to our mailing list by emailing ``listserv@listserv.ucsf.edu``
    with the message *subscribe plastidinfo firstname lastname* and an empty
    subject line

  - `Test dataset <https://www.dropbox.com/s/h17go7tnas4hpby/plastid_test_data.tar.bz2?dl=0>`_,
    for development or validation of installations
