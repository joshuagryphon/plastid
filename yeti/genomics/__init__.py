#!/usr/bin/env python
"""Defines object types useful for genomic analyses, and supplies readers
that construct these objects from data in various file formats.

    =============================================  ===========================================
    **Sub-packages**
    ------------------------------------------------------------------------------------------
    |readers|                                      Parsers for various file formats
    |scriptlib|                                    Tools for writing command-line scripts built
                                                   with ``yeti``
    ---------------------------------------------  -------------------------------------------

    **Modules**
    ------------------------------------------------------------------------------------------
    :py:mod:`~yeti.genomics.genome_array`   Randomly-accessible dictionary-like objects
                                                   that map sequencing reads and counts to
                                                   genomic positions
    :py:mod:`~yeti.genomics.genome_hash`    Randomly-accessible dictionary-like objects
                                                   that hash genomic features by position,
                                                   making looking and comparisons more efficient
    :py:mod:`~yeti.genomics.roitools`        Objects that model genomic features, even
                                                   if they are discontinuous on the chromosome
                                                   (e.g. transcripts, gapped alignments, et c).
    :py:mod:`~yeti.genomics.seqtools`       Functions for manipulating nucleic acid
                                                   sequences
    :py:mod:`~yeti.genomics.splicing`       Functions for maniuplating splice junctions
    =============================================  ===========================================
"""