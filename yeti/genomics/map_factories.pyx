# cython: embedsignature=True
"""Cython implementations of :term:`mapping rules <mapping rule>`.

See also
--------
Module documentation for :mod:`~yeti.genomics.genome_array`
    For a detailed discussion of :term:`mapping rules <mapping rule>`.
"""

import warnings
import numpy as np
cimport numpy as np
cimport c_roitools

from pysam.calignmentfile cimport AlignedSegment
from c_roitools cimport forward_strand, reverse_strand, unstranded
from yeti.util.services.exceptions import DataWarning

#from pysam.calignedsegment cimport AlignedSegment
#cdef pysam.calignedsegment.AlignedSegment AlignedSegment
cdef c_roitools.GenomicSegment GenomicSegment

INT    = np.int
FLOAT  = np.float
DOUBLE = np.double
LONG   = np.long

ctypedef np.int_t    INT_t
ctypedef np.float_t  FLOAT_t
ctypedef np.double_t DOUBLE_t
ctypedef np.long_t   LONG_t

#===============================================================================
# Factories for mapping functions for BAMGenomeArray or other structures
# Each factory returns a function that takes a list of pysam.AlignedSegments
# and a GenomicSegment, and returns a list of pysam.AlignedSegments that mapped
# to that interval under the mapping rules specified in that function, and
# a corresponding vector of read counts mapping to each nucleotide position
# in the GenomicSegment
#===============================================================================

cdef class cCenterMapFactory(object):
    """Read mapping tool for :meth:`BAMGenomeArray.set_mapping`.
    A user-specified number of bases is removed from each side of each
    read alignment, and the `N` remaining bases are each apportioned `1/N`
    of the read count, so that the entire read is counted once.

    Parameters
    ----------
    nibble : int >= 0, optional
        Number of bases to remove from each side of read (Default: `0`)
    """

    cdef unsigned int nibble

    def __cinit__(self, unsigned int nibble = 0):
        """Read mapping tool for :meth:`BAMGenomeArray.set_mapping`.
        A user-specified number of bases is removed from each side of each
        read alignment, and the `N` remaining bases are each apportioned `1/N`
        of the read count, so that the entire read is counted once.

        Parameters
        ----------
        nibble : int >= 0, optional
            Number of bases to remove from each side of read (Default: `0`)
        """
        if nibble < 0:
            raise ValueError("CenterMapFactory: `nibble` must be >= 0")

        self.nibble = nibble

    def __call__(self, list reads not None, c_roitools.GenomicSegment seg not None):
        cdef bint do_warn = 0

        cdef long seg_start = seg.start
        cdef long seg_end   = seg.end
        cdef long seg_len   = seg_end - seg_start

        # count array we will return
        cdef np.ndarray[DOUBLE_t,ndim=1] count_array = np.zeros(seg_len,dtype=DOUBLE)
        cdef double [:] count_view = count_array

        # replace by cpp vector?
        cdef list reads_out = []
        cdef list read_positions
 
        cdef:
            AlignedSegment read 
            int coord, map_length
            unsigned int read_length, i
            DOUBLE_t val
       
        for read in reads:
            read = <AlignedSegment>read
            read_positions = <list>read.positions
            read_length = len(read_positions)
            map_length = read_length - 2*self.nibble
            if map_length < 0:
                do_warn = 1
                continue
            elif map_length > 0:
                val = 1.0 / <DOUBLE_t>map_length
                for i in range(self.nibble,read_length - self.nibble):
                    coord = <long>(read_positions[i]) - seg_start
                    if coord >= 0 and coord < seg_len:
                        count_view[coord] += val

                reads_out.append(read)

        if do_warn == 1:
            warnings.warn("File contains read alignments shorter (%s nt) than `2*'nibble'` value of %s nt. Ignoring these." % (read_length,2*self.nibble),
                  DataWarning)
        
        return reads_out, count_array






 
cdef class cFivePrimeMapFactory:

    cdef int offset

    def __cinit__(self, int offset = 0):
        self.offset = offset

    def __call__(self, object reads, c_roitools.GenomicSegment seg):
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
     
 
 
cdef class f2:

    cdef int offset

    def __cinit__(self, int offset = 0):
        self.offset = offset

    def __call__(self, object reads, c_roitools.GenomicSegment seg):
        cdef:
            long seg_start = seg.start
            long seg_end   = seg.end
            long seg_len   = seg_end - seg_start
            np.ndarray[np.int_t,ndim=1] count_array = np.zeros(seg_len,dtype=int)
            list reads_out = []
            list read_positions
            AlignedSegment read
            long p_site
            int read_length
            int do_warn = 0
            int read_offset = self.offset

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
 
