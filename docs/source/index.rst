:data:`plastid` |version| welcome!
==================================


Introduction
------------

:data:`plastid` is a `Python`_ library for genomic analysis -- in particular,
:term:`high-throughput sequencing` data -- with an emphasis on simplicity for
users. It was written by Joshua Dunn in `Jonathan Weissman's lab <http://weissmanlab.ucsf.edu>`_
at `UCSF <http://ucsf.edu>`_,  initially for analysis of
:term:`ribosome profiling` and :term:`RNA-seq` data. Versions of it have been used
in several publications.

:data:`plastid`'s intended audience includes computational and traditional biologists,
software developers, and even those who are new to sequencing analysis. It is
released under the :doc:`BSD 3-Clause license <license>`.

This package provides:

  #. A set of :mod:`scripts <plastid.bin>` that implement common sequencing
     analyses

  #. A :ref:`set of classes <tour-data-structures>` that create a simple,
     intuitive interfaces to genomic :term:`features <feature>`,
     :term:`read alignments`, and quantitative data. These objects readily
     interace with existing scientific tools, like the `SciPy`_ stack.

  #. :mod:`Script writing tools <plastid.util.scriptlib>` that make it easy to
     use the :mod:`file parsers <plastid.readers>` and
     :ref:`objects <tour-data-structures>` implemented in :data:`plastid`.

  #. Extensive documentation, both in this web site and in source code


Motivation
----------
At the time of :data:`plastid`'s inception, few Python packages for genomic analysis
existed. Now, there are quite a few useful, :doc:`similar projects </related>`
-- for example, `HTSeq`_ and `metaseq`_ -- whose design goals are
similar to :data:`plastid`'s. Users are encouraged to look at those packages,
and use the package that best suits their experimental needs.

Among these, :data:`plastid` offers a number of unique features, including:

  - Classes that can represent or operate on *discontinuous* genomic features,
    like transcripts or gapped alignments, in addition to their continuous components
    (e.g. exons). See :ref:`tour-segment-chain` for more info.

  - Tools for positional analysis of sequencing data, such as the use
    of :term:`mapping rules <mapping rule>` to assign
    :term:`read alignments` to genomic positions that are *functions*
    of the alignments themselves 
    (e.g. :term:`ribosome-protected footprints <footprint>` to their P-sites),
    rather than just their 5' or 3' ends

  - Tools for easily :doc:`masking genomic positions from analysis </examples/using_masks>`,
    without needing to subtract positions from :term:`features <feature>` that
    cover those positions. This enables feature annotations to stay intact
    for analyses that don't require excluding such positions


Where to go next
----------------

Documentation & help are written for users at multiple levels of experience:

  - **Those new to sequencing and/or bioinformatics**, and those who are
    :term:`ribosome profiling` should start with :doc:`quickstart`, and then
    continue to the :doc:`tour` and :doc:`examples`. The
    :mod:`description of command-line scripts <plastid.bin>` may also be helpful.

  - **Advanced users** might be more interested in a quick :doc:`tour` of
    the primary data structures and the
    :doc:`module documentation <generated/plastid>`.

--------------------------------------------------------------------------------

Site map
========
  - :ref:`genindex`
  - :ref:`modindex`


.. _toc-getting-started:
 
.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Getting started
    
   quickstart
   tour
   installation
   List of command-line scripts </scriptlist>


.. _toc-user-manual:
 
.. toctree::
   :hidden:
   :maxdepth: 5
   :caption: User manual
    
   examples
   Module documentation <generated/plastid>
   FAQ
   glossary
   zreferences


.. _toc-developer:

.. toctree::
   :hidden:
   :includehidden:
   :maxdepth: 10
   :caption: Developer info

   Entrypoints </devinfo/entrypoints>
   Contributing </devinfo/contributing>

.. _toc-other-info:
     
.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Other information
    
   citing
   License <license>
   CHANGES
   related
   contact
