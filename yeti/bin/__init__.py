#!/usr/bin/env python
"""Command-line scripts that implement common sequencing workflows, provided
a genome :term:`annotation` file and/or a file containing :term:`read alignments`

    ======================   ==========================================================
    **Sequencing analysis**                   
    -----------------------------------------------------------------------------------
    |counts_in_region|       Count the number of :term:`read alignment` s covering
                             arbitrary regions of interest in the genome, and calculate
                             read densities (in reads per nucleotide and in :term:`RPKM`)
                             over these regions
    
    |cs|                     Count the number of :term:`read alignments<alignment>`
                             and calculate read densities (in :term:`RPKM`)
                             specifically for genes and sub-regions (5' UTR,
                             CDS, 3' UTR)

    |get_count_vectors|      Fetch vectors of :term:`counts` at each nucleotide position
                             in one or more regions of interest, saving each vector
                             as its own line-delimited text file
                             
    |make_wiggle|            Create `wiggle`_ or `bedGraph`_ files from alignment files
                             after applying a read :term:`mapping rule` (e.g.
                             to map :term:`ribosome-protected footprints <footprint>`
                             at their :term:`P-sites <P-site offset>`), for
                             visualization in a :term:`genome browser`
                             
    |metagene|               Compute a :term:`metagene` profile of :term:`read alignments`
                             or :term:`counts` over one or more regions of interest,
                             optionally applying a :term:`mapping rule`
                             
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