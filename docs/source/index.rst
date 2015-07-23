:data:`yeti` |version| welcome!
===============================

Introduction
------------

:py:data:`yeti` is a `Python`_ library for genomic analysis -- in particular, :term:`high-throughput sequencing` data -- with an emphasis on simplicity for users.

Its intended audience includes computational and traditional biologists, software developers, and even those who are new to sequencing analysis. It is released under the :doc:`BSD 3-Clause license <license>`.

This package provides:

  #. A set of :mod:`scripts <yeti.bin>` that implement common sequencing analyses

  #. A :ref:`set of objects <tour-data-structures>` that create a simple, intuitive interfaces to genomic :term:`features <feature>`, :term:`read alignments`, and quantitative data. These objects readily interace with existing scientific tools, like the `SciPy`_ stack.

  #. :mod:`Script writing tools <yeti.util.scriptlib>` that make it easy to use the :mod:`file parsers <yeti.readers>` and :ref:`objects <tour-data-structures>` implemented in :data:`yeti`.

s
Where to go next
----------------

Documentation & help are written for users at multiple levels of experience.

  - **Those new to sequencing and/or bioinformatics**, and those who are :term:`ribosome profiling` should start with :doc:`quickstart`, and then continue to the :doc:`tour` and :doc:`examples`. The :mod:`description of command-line scripts <yeti.bin>` may also be helpful.

  - **Advanced users** might be more interested in a quick :doc:`tour` of the primary data structures and the :doc:`module documentation <generated/yeti>`.

--------------------------------------------------------------------------------

Site map
========
  - :ref:`genindex`
  - :ref:`modindex`


.. _toc-getting-started:
 
.. toctree::
   :maxdepth: 2
   :caption: Getting started
    
   quickstart
   tour
   installation


.. _toc-user-manual:
 
.. toctree::
   :maxdepth: 3
   :caption: User manual
    
   examples
   FAQ
   glossary
   zreferences


.. _toc-technical-documentation:
   
.. toctree::
   :includehidden:
   :maxdepth: 4
   :caption: Technical documentation
   
   Module documentation <generated/yeti>


.. _toc-other-info:
     
.. toctree::
   :maxdepth: 1
   :caption: Other information
    
   citing
   License <license>
   contributing
   CHANGES
   related
   contact
