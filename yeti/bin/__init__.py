#!/usr/bin/env python
"""Command-line scripts that implement common sequencing workflows

    =========================   =============================================================================
    **Analysis of sequencing or quantitative data**                 
    ---------------------------------------------------------------------------------------------------------
    |counts_in_region|           Count the number of :term:`read alignments <alignment>` covering
                                 arbitrary regions of interest in the genome, and calculate
                                 read densities (in reads per nucleotide and in :term:`RPKM`)
                                 over these regions
    
    |cs|                         Count the number of :term:`read alignments<alignment>`
                                 and calculate read densities (in :term:`RPKM`)
                                 specifically for genes and sub-regions (5' UTR,
                                 CDS, 3' UTR)

    |get_count_vectors|          Fetch vectors of :term:`counts` or other quantitative data
                                 at each nucleotide position
                                 in one or more regions of interest, saving each vector
                                 as its own line-delimited text file
                             
    |make_wiggle|                Create `wiggle`_ or `bedGraph`_ files from alignment files
                                 after applying a read :term:`mapping rule` (e.g.
                                 to map :term:`ribosome-protected footprints <footprint>`
                                 at their :term:`P-sites <P-site offset>`), for
                                 visualization in a :term:`genome browser`
                             
    |metagene|                   Compute a :term:`metagene` profile of :term:`read alignments`,
                                 :term:`counts`, or quantitative data over one or more regions of interest,
                                 optionally applying a :term:`mapping rule`
                             
    |phase_by_size|              Estimate :term:`sub-codon phasing` in
                                 :term:`ribosome profiling` data
    
    |psite|                      Estimate position of ribosomal P-site within
                                 :term:`ribosome profiling` :term:`read alignments`
                                 as a function of read length
    -------------------------   -----------------------------------------------------------------------------
    **Generating or modifying genome annotations**    
    ---------------------------------------------------------------------------------------------------------
    |crossmap|                   Generate a :term:`mask file` annotating regions of the genome
                                 that fail to produce
                                 uniquely mapping reads under various alignment
                                 criteria, so that they may be excluded from future analyses
                             
    |gff_parent_types|           Examine parent-child relationships of features
                                 in a `GFF3`_ file
                             
    |reformat_transcripts|       Convert transcripts between `BED`_, `BigBed`_,
                                 `GTF2`_, `GFF3`_, and `PSL`_ formats
                             
    |findjuncs|                  Find all unique splice junctions in one or more
                                 transcript annotations, and optionally export these
                                 in `Tophat`_'s ``.juncs`` format
                             
    |slidejuncs|                 Compare splice junctions discovered in data to those
                                 in an annotation of known splice junctions, and,
                                 if possible with equal sequence support, slide
                                 discovered junctions to compatible known junctions.
    =========================   =============================================================================
    
"""