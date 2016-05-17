#!/usr/bin/env python
"""This module contains a parser for `bowtie`_'s legacy output format.
Functions in this module are are seldom used on their own, and rather are
accessed by |GenomeArray| or |SparseGenomeArray| when importing from
`bowtie`_ files.

See also
--------
|GenomeArray| and |SparseGenomeArray|
    Array-like objects that store and index quantitative data over genomes

http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-output
    Detailed description of `bowtie`_ output format

:py:mod:`plastid.test.unit.genomics.test_genome_array`
    for integrative tests of these functions
"""
__author__ = "joshua"
__date__ = "2011-03-18"

from plastid.genomics.roitools import SegmentChain, GenomicSegment
from plastid.util.io.filters  import AbstractReader


#===============================================================================
# INDEX: Readers for various bowtie1-like alignment file formats
#===============================================================================


class BowtieReader(AbstractReader):
    """Read alignments from `bowtie`_ files line-by-line into |SegmentChains|.
    The following attributes are defined and stored in the `attr` dict of
    each returned |SegmentChain|
    
    `seq_as_aligned`
        the sequence in the direction it aligns, NOT necessarily
        the read in the direction it was sequenced
    
    `qualstr_phred`
        a quality string, phred encoded
    
    `total_alignments`
        the number of total alignments found
    
    See description of `bowtie`_ legacy format at
    http://bowtie-bio.sourceforge.net/manual.shtml
    
    Parameters
    ----------
    stream : file-like
        Stream of alignments in `bowtie`_'s legacy output format
        
    
    Yields
    -------
    |SegmentChain|
        A read alignment
    """
    def filter(self,line):
        """Parse a read alignment as |SegmentChain| from a line of `bowtie`_ output"""
        items = line.strip("\n").split("\t")
        read_name      = items[0]
        strand         = items[1]
        ref_seq        = items[2]
        coord          = int(items[3])
        attr = { 'seq_as_aligned' : items[4],
                 'qualstr'        : items[5],
                 'mismatch_str'   : items[7],
                 'type'           : "alignment",
                 'ID'             : read_name,
               }
        
        iv = GenomicSegment(ref_seq,coord,coord+len(attr['seq_as_aligned']),strand)
        feature = SegmentChain(iv,**attr)
        return feature
