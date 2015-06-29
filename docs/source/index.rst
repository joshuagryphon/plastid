:data:`yeti` v0.5 welcome!
==========================

Introduction
------------

:py:data:`yeti` is a lightweight Python library for analysis of 
:term:`high-throughput sequencing` data. Its intended audience
includes computational as well as traditional biologists, even those
who may be new to sequencing analysis.

To this end, :data:`yeti` seeks to flatten the learning curve required to
analyze genomics data interactively. Our design goals are to:

  - simplify access to data, regardless of its underlying format. To do so,
    :data:`yeti` provides a unified and intuitive set of interfaces to:

      - annotation data in `BED`_, `BigBED`_, `GTF2`_, `GFF3`_, or `PSL`_ format

      - quantitative data in `Wiggle`_ or `bedGraph`_ format

      - read alignments in `BAM`_ (via `Pysam`_) or `bowtie's native format <bowtie>`_

  - easily integrate into the Python ecosystem, especially the
    `SciPy stack <http://www.scipy.org/stackspec.html>`_

  - provide a set of tools to data nucleotide-by-nucleotide over a region of
    interest, instead of just, for example, counting the number of bulk read
    alignments that cross a region

  - do all of this as simply as possible



Package contents
----------------

:data:`yeti` includes:

  - a number of :ref:`scripts <scripts>` that implement common sequencing analyses

  - a `code library <generated/yeti>`_ of :doc:`data structures <overview>` for
    interactive analysis and `scripting components </generated/yeti/util/scriptlib>`_
    to simplify writing new scripts

  - a comprehensive test suite, to keep bugs to a minimum



Where to go next
----------------

**Those new to sequencing**, and those who are :term:`ribosome profiling`
should start with :doc:`quickstart`, and then continue to the :ref:`cookbook`
and/or :ref:`scripts`.

**Advanced users** might be more interested in a quick :ref:`oveview <overview>`, 
and the `technical documentation <generated/yeti>`_


   
.. Indices and tables
.. ------------------

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


 .. toctree::
    :maxdepth: 2
