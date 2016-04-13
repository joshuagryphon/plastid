:data:`plastid` |version| welcome!
==================================


Introduction
------------

:data:`plastid` is a `Python`_ library for genomics and sequencing that seeks to 
flatten the learning curve for both bench and computational biologists. It
:ref:`intuitive tools for exploratory data analysis <tour-data-structures>` (EDA),
as well as :mod:`command-line scripts <plastid.bin>` that implement standard
analyses.

For EDA, plastid defines consistent interfaces for the many similar
file types that exist in genomics, allowing users to focus on biology, without
first having to understand the quirks of file formats. The interfaces themselves
use familiar biological idioms, such as :ref:`spliced transcripts <tour-segment-chain>`.

In addition, plastid's data models interface directly with the rich computing environment
provided by the `SciPy`_ stack, enabling development of sophisticated analyses with
substantially reduced effort.

:data:`plastid` was written by Joshua Dunn in `Jonathan Weissman's lab <http://weissmanlab.ucsf.edu>`_
at `UCSF <http://ucsf.edu>`_. Versions of it have been used in several publications.



How it's unique
---------------

:data:`plastid` offers a number of unique features, including:

 - handling discontinuous genomic features, such as transcripts or gapped alignments,
   in a single object, called a :ref:`SegmentChain <tour-segment-chain>`

 - providing a suite of tools for positional analysis of sequencing data, such as

    - :term:`mapping rules <mapping rule>` that can exploit arbitrary properties of a 
      :term:`read alignment <read alignments>` to determine the genomic position(s)
      to which it should map 
      (e.g. map :term:`ribosome-protected footprints <footprint>` to their P-sites,
      rather than just their 5' or 3' ends)

    - tools for  :doc:`masking genomic positions from analysis </examples/using_masks>`,
      without needing to subtract positions from :term:`features <feature>` that
      cover those positions. This enables feature annotations to stay intact
      for analyses that don't require excluding such positions

 - accepting simple plugin functions (like :term:`mapping rules <mapping rule>`) that
   globally modify its behavior




Where to go next
----------------

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

   Contributing </devinfo/contributing>
   Entrypoints </devinfo/entrypoints>

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
