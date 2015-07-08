:data:`yeti` v0.5 welcome!
==========================

Introduction
------------

:py:data:`yeti` is a `Python`_ library for genomic analysis -- in particular,
:term:`high-throughput sequencing` data -- with an emphasis on simplicity
for users.

Its intended audience includes computational biologists, traditional biologists,
software developers, and even those who may be new to sequencing analysis.

This package provides:

  #. :mod:`Command-line scripts <yeti.bin>` that implement common sequencing
     workflows
  
  #. APIs that simplify access to genomics data -- regardless of the underlying formats -- 
     and readily interface with scientific tools like the `SciPy`_ stack.
     These include a unified and intuitive set of interfaces to:

      - annotation data in `BED`_, `BigBED`_, `GTF2`_, `GFF3`_, or `PSL`_ format

      - quantitative data in `Wiggle`_ or `bedGraph`_ format

      - read alignments in `BAM`_ (via `Pysam`_) or `bowtie's native format <bowtie>`_
     
     Readers for these file formats may be found in the |readers| subpackage.
     The API documentation for the objects they create may be found in the
     |genomics| subpackage.

  #. :mod:`Script writing tools <yeti.util.scriptlib>` that make it easy to implement
     new workflows as command-line scripts.

  #. A comprehensive :mod:`test suite <yeti.test>`, to keep bugs to a minimum


Where to go next
----------------

Documentation & help are written for users at multiple levels of experience.

  * **Those new to sequencing and/or bioinformatics**, and those who are
    :term:`ribosome profiling` should start with :doc:`quickstart`, and then
    continue to the :ref:`cookbook` and/or :ref:`scripts`. The :doc:`concepts
    <concepts>` section may also be helpful.

  * **Advanced users** might be more interested in a quick :ref:`overview <overview>`
    of the relevant data structures and the `technical documentation <generated/yeti>`_.

   
Indices, tables, and links
--------------------------

  * `our github link`_
  * `pypi link`_
  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`


 .. toctree::
