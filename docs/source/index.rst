:data:`plastid` |version|
=========================


Welcome!
--------

:data:`plastid` is a `Python`_ library for genomics and sequencing. It includes
:ref:`tools for exploratory data analysis <tour-data-structures>` (EDA) 
as well as a handful of :doc:`scripts </scriptlist>` that implement common tasks.

:data:`plastid` differs from other packages in its design goals. Namely:

 - its intended audience includes both **bench and computational biologists**.
   We tried to make it easy to use, and wrote lots of :doc:`/examples`
 
 
 - It is designed for analyses in which data at each position
   **within a gene or transcript** are of interest, such as analysis of
   :term:`ribosome profiling` data. To this end, :data:`plastid`

    - uses :doc:`/concepts/mapping_rules` to **extract the biology of interest**
      from :term:`read alignments` -- e.g. in the case of 
      :term:`ribosome profiling`, a ribosomal :term:`P-site`, in :term:`DMS-seq`,
      sites of nucleotide modification, et c. --
      and turn these into quantitative data, usually
      :class:`numpy arrays <numpy.ndarray>` of counts at each nucleotide
      position in a transcript.
  
    - **encapsulates multi-segment features**, such as spliced transcripts,
      as :ref:`single objects <tour-segment-chain>`. This facilitates many common
      tasks, such as converting coordinates between genome and feature-centric
      spaces.
   
 
 - It **separates data from its representation** on disk by providing consistent
   interfaces to many of the various :ref:`file formats <file-format-table>`, 
   found in the wild.
   
   
 - It is **designed for expansion** to new or unknown assays. Frequently,
   :ref:`writing a new mapping rule <mapping-rules-roll-your-own>` is sufficient
   to enable all of :data:`plastid`'s tools to interpret data coming from a new
   assay.



:data:`plastid` was written by Joshua Dunn in
`Jonathan Weissman's lab <http://weissmanlab.ucsf.edu>`_ at
`UCSF <http://ucsf.edu>`_. Versions of it have been used in several publications
(:cite:`Dunn2013,FieldsRodriguez2015`).



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
--------
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
   Demo dataset </test_dataset>
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
