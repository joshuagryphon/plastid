#!/usr/bin/env python
"""Command-line scripts that implement common sequencing workflows.

    ======================   ==========================================================
    **Sequencing analysis**                   
    -----------------------------------------------------------------------------------
    |counts_in_region|       Compute raw counts and RPKM for arbitrary regions
                             of interest
    |cs|                     Compute raw counts and RPKM for genes, transcripts,
                             5'UTRs, CDS, and 3'UTRs, excluding positions overlapped
                             by neighboring genes
    |get_count_vectors|      Compute vectors of counts at each position in one
                             or more regions of interest
    |make_wiggle|            Create `wiggle`_ or `bedGraph`_ files from alignments
                             after applying read mapping (e.g. P-site mapping) rules
    |metagene|               Compute metagene averages of read profiles over one
                             or more regions of interest, after applying read
                             mapping rules
    |phase_by_size|          Compute sub-codon periodicity for ribosome profiling data
    |psite|                  Estimate position of ribosomal P-site within ribosome
                             profiling reads, as a function of read length
    ----------------------   ----------------------------------------------------------
    **Generating or modifying genome annotations**    
    -----------------------------------------------------------------------------------
    |crossmap|               Determine which regions of the genome fail to produce
                             uniquely mapping reads under various alignment
                             criteria, so that these regions may be excluded
                             from downstream analyses
    |gff_parent_types|       Determine parent-child relationships for features
                             in a GFF3 file
    |reformat_transcripts|   Convert transcripts between `BED`_, `BigBed`_,
                             `GTF2`_, `GFF3`_, and `PSL`_ formats
    |findjuncs|              Find all unique splice junctions in one or more
                             transcript annotations, and optionally export these
                             in `Tophat`_'s ``.juncs`` format
    |slidejuncs|             Compare splice junctions discovered in data to those
                             in an annotation of known splice junctions, and,
                             if possible with equal sequence support, slide
                             discovered junctions to compatible known junctions.
    ======================   ==========================================================
    
"""