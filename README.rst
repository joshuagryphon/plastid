Welcome to `plastid`!
=====================

For documentation, see `our home page <http://plastid.readthedocs.org/en/latest/>`_
on `ReadtheDocs.org <http://readthedocs.org>`_.

To run the tests, download the `test dataset <https://www.dropbox.com/s/h17go7tnas4hpby/plastid_test_data.tar.bz2?dl=0>`_ and unpack
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
released under the BSD 3-Clause license.

This package provides:

  #. A set of scripts that implement common sequencing analyses,

  #. A set of classes that create simple, intuitive interfaces to complex
     genomic features, read alignments, and quantitative data. These objects
     readily interface with existing scientific tools, like the SciPy stack,
     to facilitate interactive / exploratory data analysis.

  #. Script writing tools that make it easy to use the objects
     implemented in ``plastid``.

  #. Extensive documentation, both in source code and at 
     `our home page <http://plastid.readthedocs.org/en/latest/>`_
     on `ReadtheDocs.org <http://readthedocs.org>`_.


Installation
------------
``plastid`` can be installed directly from PyPI.

    1. First make sure you have numpy and pysam installed::

        pip install numpy pysam

    2. Then::
        
        pip install plastid

    3. If when running you get any warnings about numpy versions having changed,
       regenerate the c files to match your numpy headers by typing:

       pip install cython
       pip install plastid --recythonize


Links & help
------------

  - `Documentation <http://plastid.readthedocs.org>`_

  - `Our github repo <https://github.com/joshuagryphon/plastid`_

  - Subscribe to our mailing list by emailing ``listserv@listserv.ucsf.edu``
    with the message *subscribe plastidinfo firstname lastname* and
    an empty subject line.

  - Test dataset <https://www.dropbox.com/s/h17go7tnas4hpby/plastid_test_data.tar.bz2?dl=0>`_,
    for developers or validating installs
