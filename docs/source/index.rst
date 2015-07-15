:data:`yeti` |version| welcome!
===============================

Introduction
------------

:py:data:`yeti` is a `Python`_ library for genomic analysis -- in particular,
:term:`high-throughput sequencing` data -- with an emphasis on simplicity
for users.

Its intended audience includes computational biologists, traditional biologists,
software developers, and even those who may be new to sequencing analysis. It
is released under the :doc:`BSD 3-Clause license <license>`.

This package provides:

  #. A set of :mod:`scripts <yeti.bin>` that implement common sequencing
     analyses
  
  #. A set of objects that create a simple and intuitive interface
     to various types of data used in genomics. These objects
     readily interace with existing scientific tools, like the
     `SciPy`_ stack. At present, :data:`yeti` can read:

      - annotation data in `BED`_, `BigBED`_, `GTF2`_, `GFF3`_, or `PSL`_ format,
        with or without `tabix`_ compression

      - quantitative data in `Wiggle`_ or `bedGraph`_ format

      - read alignments in `BAM`_ (via `Pysam`_) or `bowtie's native format <bowtie>`_
     
     Readers for these file formats may be found in the |readers| subpackage.
     The API documentation for the objects they create may be found in the
     |genomics| subpackage.

  #. :mod:`Script writing tools <yeti.util.scriptlib>` that make it easy to implement
     new workflows as command-line scripts.


Where to go next
----------------

Documentation & help are written for users at multiple levels of experience.

  * **Those new to sequencing and/or bioinformatics**, and those who are
    :term:`ribosome profiling` should start with :doc:`quickstart`, and then
    continue to the :doc:`examples` and/or
    :mod:`description of command-line scripts <yeti.bin>`. :doc:`concepts`
    section may also be helpful.

  * **Advanced users** might be more interested in a quick :doc:`tour <tour>`
    of the primary data structures and the :doc:`module documentation <generated/yeti>`.

   
Indices and links
-----------------

  - :doc:`Module documentation <generated/yeti>`
  - :ref:`Module index <modindex>`
  - :ref:`Alphabetical index of functions and classes <genindex>`
  - `Package download <pypi link>`_ from `PyPI <http://pypi.python.org>`_
  - `Source code for developers or contributors, on github <our github link>`_

.. .. toctree::
