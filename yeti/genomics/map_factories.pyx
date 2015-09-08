import warnings
import numpy as np
cimport numpy as np

cimport c_roitools
from c_roitools import forward_strand, reverse_strand, unstranded

cimport pysam.calignedsegment
from yeti.util.services.exceptions import DataWarning

from pysam.calignedsegment cimport AlignedSegment

#cdef pysam.calignedsegment.AlignedSegment AlignedSegment
cdef c_roitools.GenomicSegment GenomicSegment

#===============================================================================
# Factories for mapping functions for BAMGenomeArray or other structures
# Each factory returns a function that takes a list of pysam.AlignedSegments
# and a GenomicSegment, and returns a list of pysam.AlignedSegments that mapped
# to that interval under the mapping rules specified in that function, and
# a corresponding vector of read counts mapping to each nucleotide position
# in the GenomicSegment
#===============================================================================

cdef class cCenterMapFactory:

    cdef int nibble

    def __cinit__(self, int nibble = 0):
        self.nibble = nibble

    def __call__(self, list reads, c_roitools.GenomicSegment seg):
        cdef long seg_start = seg.start
        cdef long seg_end   = seg.end
        cdef long seg_len   = seg_end - seg_start
        cdef frozenset seg_positions = frozenset(range(seg_len))
        cdef frozenset overlap
        cdef list reads_out = []
        cdef list read_positions
        cdef list overlap_array_coordinates
        cdef np.ndarray[np.double_t,ndim=1] count_array = np.zeros(len(seg),dtype=float)
        cdef AlignedSegment read
        cdef int read_length
        cdef double val
        cdef int do_warn = 0

        for read in reads:
            read_positions = read.positions
            read_length = len(read_positions)
            if read_length <= 2*self.nibble:
                do_warn = 1
                continue
            
            if self.nibble != 0:
                read_positions = read_positions[self.nibble:-self.nibble]
 
            overlap = frozenset(read_positions) & seg_positions
            if len(overlap) > 0:
                overlap_array_coordinates = [X-seg_start for X in overlap]
                reads_out.append(read)
                val = 1.0 / read_length
                count_array[overlap_array_coordinates] += val

        if do_warn == 1:
            warnings.warn("File contains read alignments shorter (%s nt) than `2*'nibble'` value of %s nt. Ignoring these." % (read_length,2*self.nibble),
                  DataWarning)
        
        return reads_out, count_array
 
 
cdef class cFivePrimeMapFactory:

    cdef int offset

    def __cinit__(self, int offset = 0):
        self.offset = offset

    def __call__(self, list reads, c_roitools.GenomicSegment seg):
        cdef long seg_start = seg.start
        cdef long seg_end   = seg.end
        cdef long seg_len   = seg_end - seg_start
        cdef np.ndarray[np.int_t,ndim=1] count_array = np.zeros(seg_len,dtype=int)
        cdef list reads_out = []

        cdef list read_positions
        cdef AlignedSegment read
        cdef long p_site
        cdef int read_length
        cdef int do_warn = 0
        
        cdef int read_offset = self.offset
        if seg.strand == reverse_strand:
             read_offset = -read_offset - 1 

        for read in reads:
            read_positions = read.positions
            read_length = len(read_positions)
            if self.offset >= read_length:
                do_warn = 1
                continue
             
            p_site = read_positions[read_offset]
            if p_site >= seg_start and p_site < seg_end:
                reads_out.append(read)
                count_array[p_site - seg_start] += 1

        if do_warn == 1:
            warnings.warn("File contains read alignments shorter (%s nt) than offset (%s nt). Ignoring." % (read_length,self.offset),
                          DataWarning)

        return reads_out, count_array
     
 
 

