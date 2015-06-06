#!/usr/bin/env python
"""Contains utility functions for reading and writing alignment information
from bowtie's proprietary legacy output format, as well as functions that 
map read alignments to specific positions. Functions defined in this class
are are seldom used on their own, and rather are accessed internally by |GenomeArray|
when importing from bowtie or TagAlign files.

See also
--------
http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-output
    for detailed description of bowtie output format

:py:mod:`yeti.test.unit.genomics.test_genome_array`
    for integrative tests of theese functions
"""
__author__ = "joshua"
__date__ = "2011-03-18"

from yeti.genomics.roitools import SegmentChain, GenomicSegment
from yeti.util.io.filters  import AbstractReader


#===============================================================================
# INDEX: Readers for various bowtie1-like alignment file formats
#===============================================================================


class BowtieReader(AbstractReader):
    """Reads bowtie files into |SegmentChain| objects corresponding to alignments.
    Attributes defined for the feature are:
    
    `seq_as_aligned`
        the sequence in the direction it aligns, NOT necessarily
        the read in the direction it was sequenced
    
    `qualstr_phred`
        a quality string, phred encoded
    
    `total_alignments`
        the number of total alignments found
    
    See description of format at http://bowtie-bio.sourceforge.net/manual.shtml
    
    Returns
    -------
    |SegmentChain|
    """
    def filter(self,line):
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


class TagalignReader(AbstractReader):
    """Reads alignments from Nick Ingolia's TagAlign utility to |SegmentChain|
    Attributes defined for this feature are:

    `seq_as_aligned`
        the sequence in the direction it aligns, NOT necessarily
        the read in the direction it was sequenced
    
    `qualstr_phred`
        a quality string, phred encoded
    
    `total_alignments`
        the number of total alignments found

    Returns
    -------
    |SegmentChain|
    """
    def filter(self,line):
        items = line.strip("\n").split("\t")
        read_name      = items[0]
        strand         = items[8]
        chrom          = items[6]
        attr = { 'seq_as_aligned' : items[1],
                 'qualstr'        : items[2],
                 'total_alignments' : int(items[3]),
                 'type'           : "alignment",
                 'ID'             : read_name,
               }
        fiveprime_align  = int(items[7])
        max_align_length = int(items[5])
        if strand == "+":
            start = fiveprime_align
            end   = fiveprime_align + max_align_length
        elif strand == "-":
            end   = fiveprime_align + 1
            start = end - max_align_length
        iv = GenomicSegment(chrom,start,end,strand)
        feature = SegmentChain(iv,**attr)
        return feature
