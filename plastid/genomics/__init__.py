#!/usr/bin/env python
"""This package contains several object types useful for genomic analyses.

Package overview
================

    =============================================  ==================================================================
    **Submodule**                                   **Description**
    ---------------------------------------------  ------------------------------------------------------------------
    :py:mod:`~plastid.genomics.genome_array`         Array-like objects indexed by chromosome, position, and strand
                                                     
    :py:mod:`~plastid.genomics.genome_hash`          Dictionary-like objects associate features
                                                     with genomic positions, allowing fast query of which
                                                     features overlap a region or nucleotide

    :py:mod:`~plastid.genomics.map_factories`        :term:`Mapping rules <mapping rule>` for
                                                     :class:`~plastid.genomics.genome_array.BAMGenomeArray`
                                                   
    :py:mod:`~plastid.genomics.roitools`             Objects that represent genomic :term:`features <feature>`,
                                                     such as genes or transcripts
    
    :py:mod:`~plastid.genomics.seqtools`             Functions for manipulating nucleic acid
                                                     sequences
                                            
    :py:mod:`~plastid.genomics.splicing`             Functions for manipulating splice junctions
    =============================================  ==================================================================
"""
