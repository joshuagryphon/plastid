#!/usr/bin/env python
"""This package contains several object types useful for genomic analyses.

Package overview
================

    =====================================  =====================================================
    **Submodules**
    --------------------------------------------------------------------------------------------
    :py:mod:`~yeti.genomics.genome_array`   Dictionary-like objects that associate read alignments
                                            and/or quantitative data with genomic positions
                                                   
    :py:mod:`~yeti.genomics.genome_hash`    Dictionary-like objects associate genomic features
                                            with genomic positions, allowing fast query of which
                                            features overlap a region or nucleotide
                                                   
    :py:mod:`~yeti.genomics.roitools`       Objects that represent genomic  :term:`features <feature>`,
                                            such as genes or transcripts
    
    :py:mod:`~yeti.genomics.seqtools`       Functions for manipulating nucleic acid
                                            sequences
                                            
    :py:mod:`~yeti.genomics.splicing`       Functions for manipulating splice junctions
    =====================================  =====================================================
"""