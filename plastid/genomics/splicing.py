#!/usr/bin/env python
"""Utility functions for systematically naming splice junctions and retrieving their sequences
"""
import re

#===============================================================================
# INDEX: Utility functions for naming splice junctions and retrieving sequences
#===============================================================================

junction_pat = re.compile(r"([^@]+)@([0-9]+)\^([0-9]+)\(([+-])\)f([0-9]+)")

def get_junction_tuple(ivc):
    """Convert an |SegmentChain| representing a splice junction to a tuple
    
    Parameters
    ----------
    ivc : |SegmentChain|
         A two-exon fragment representing a splice junction
        
    Returns
    -------
    tuple
        `(chromosome name, half-open end of fiveprime exon, first position of threeprime exon, strand)`
    """
    return (ivc.spanning_segment.chrom,ivc[0].end,ivc[1].start,ivc.spanning_segment.strand)
