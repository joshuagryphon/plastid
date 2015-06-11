#!/usr/bin/env python
"""Defines object types useful for genomic analyses, and supplies readers
that construct these objects from data in various file formats.

    =====================================  =====================================================
    **Sub-packages**
    --------------------------------------------------------------------------------------------
    |readers|                                      Parsers for various file formats
    |scriptlib|                                    Tools for writing command-line scripts built
                                                   with ``yeti``
    -------------------------------------  -----------------------------------------------------
    **Modules**
    --------------------------------------------------------------------------------------------
    :py:mod:`~yeti.genomics.genome_array`   Randomly-accessible dictionary-like objects
                                            that map :term:`read alignments` and quantitative
                                            data to genomic positions
                                                   
    :py:mod:`~yeti.genomics.genome_hash`    Randomly-accessible dictionary-like objects
                                            that fetch genomic features that overlap
                                            or reside within genomic regions
                                                   
    :py:mod:`~yeti.genomics.roitools`       Objects that represent genomic 
                                            term:`features <feature>`
    
    :py:mod:`~yeti.genomics.seqtools`       Functions for manipulating nucleic acid
                                            
                                            sequences
    :py:mod:`~yeti.genomics.splicing`       Functions for manipulating splice junctions
    =====================================  =====================================================
"""