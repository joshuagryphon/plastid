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
            raise ValueError("CenterMapFactory: `nibble` must be >= 0. Got %s." % nibble)

        self.nibble = nibble

    def __call__(self, list reads not None, c_roitools.GenomicSegment seg not None):
        """Returns reads covering a region, and a count vector mapping reads
        to specific positions in the region. `self.nibble` bases are trimmed from each
        side of the read, and each of the `N` remaining alignment positions
        are incremented by `1/N`
        
        Parameters
        ----------
        seg : |GenomicSegment|
            Region of interest
        
        reads : list of :py:class:`pysam.AlignedSegment`
            Reads to map
            
        Returns
        -------
        :py:class:`numpy.ndarray`
            Vector of counts at each position in `seg`
        """
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
            read_positions = <list>(<AlignedSegment>read).positions
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
            warnings.warn("Data contains read alignments shorter (%s nt) than `2*'nibble'` value of %s nt. Ignoring these." % (read_length,2*self.nibble),
                  DataWarning)
        
        return reads_out, count_array


cdef class cFivePrimeMapFactory:
    """Fiveprime mapping factory for :py:meth:`BAMGenomeArray.set_mapping`.
    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
     
    Parameters
    ----------
    offset : int >= 0, optional
        Offset from 5' end of read, in direction of 3' end, at which 
        reads should be counted (Default: `0`)
    """
    cdef int offset

    def __cinit__(self, int offset = 0):
        """Fiveprime mapping factory for :py:meth:`BAMGenomeArray.set_mapping`.
        Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
         
        Parameters
        ----------
        offset : int >= 0, optional
            Offset from 5' end of read, in direction of 3' end, at which 
            reads should be counted (Default: `0`)
        """
        if offset < 0:
            raise ValueError("FivePrimeMapFactory: `offset` must be <= 0. Got %s." % offset)
        self.offset = offset

    def __call__(self, list reads not None, c_roitools.GenomicSegment seg not None):
        """Returns reads covering a region, and a count vector mapping reads
        to specific positions in the region, mapping reads at `self.offset`
        from the fiveprime end of each read.
 
        Parameters
        ----------
        reads : list of :py:class:`pysam.AlignedSegment`
            Reads to map
            
        seg : |GenomicSegment|
            Region of interest

            
        Returns
        -------
        list
            List of :py:class:`pysam.AlignedSegment` that map to region
            
        :py:class:`numpy.ndarray`
            Vector of counts at each position in `seg`
        """
        cdef:
            long seg_start = seg.start
            long seg_end   = seg.end
            long seg_len   = seg_end - seg_start
            np.ndarray[LONG_t,ndim=1] count_array = np.zeros(seg_len,dtype=LONG)
            long [:] count_view = count_array

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
            read_positions = <list>read.positions
            read_length = len(read_positions)
            if self.offset >= read_length:
                do_warn = 1
                continue
             
            p_site = read_positions[read_offset]
            if p_site >= seg_start and p_site < seg_end:
                reads_out.append(read)
                count_view[p_site - seg_start] += 1

        if do_warn == 1:
            warnings.warn("Data contains read alignments shorter (%s nt) than offset (%s nt). Ignoring." % (read_length,self.offset),
                          DataWarning)

        return reads_out, count_array
     
cdef class cThreePrimeMapFactory:
    """Threeprime mapping factory for :py:meth:`BAMGenomeArray.set_mapping`.
    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
     
    Parameters
    ----------
    offset : int >= 0, optional
        Offset from 3' end of read, in direction of 5' end, at which 
        reads should be counted (Default: `0`)
    """
    cdef int offset

    def __cinit__(self, int offset = 0):
        """Threeprime mapping factory for :py:meth:`BAMGenomeArray.set_mapping`.
        Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
         
        Parameters
        ----------
        offset : int >= 0, optional
            Offset from 3' end of read, in direction of 5' end, at which 
            reads should be counted (Default: `0`)
        """
        if offset < 0:
            raise ValueError("ThreePrimeMapFactory: `offset` must be <= 0. Got %s." % offset)
        self.offset = offset

    def __call__(self, list reads not None, c_roitools.GenomicSegment seg not None):
        """Returns reads covering a region, and a count vector mapping reads
        to specific positions in the region, mapping reads at `self.offset`
        from the threeprime end of each read.
 
        Parameters
        ----------
        reads : list of :py:class:`pysam.AlignedSegment`
            Reads to map
            
        seg : |GenomicSegment|
            Region of interest

            
        Returns
        -------
        list
            List of :py:class:`pysam.AlignedSegment` that map to region
            
        :py:class:`numpy.ndarray`
            Vector of counts at each position in `seg`
        """
        cdef:
            long seg_start = seg.start
            long seg_end   = seg.end
            long seg_len   = seg_end - seg_start
            np.ndarray[LONG_t,ndim=1] count_array = np.zeros(seg_len,dtype=LONG)
            long [:] count_view = count_array

            list reads_out = []
            list read_positions
            AlignedSegment read
            long p_site
            int read_length
            int do_warn = 0
            int read_offset = self.offset

        if seg.strand != reverse_strand:
             read_offset = -read_offset - 1 

        for read in reads:
            read_positions = <list>read.positions
            read_length = len(read_positions)
            if self.offset >= read_length:
                do_warn = 1
                continue
             
            p_site = read_positions[read_offset]
            if p_site >= seg_start and p_site < seg_end:
                reads_out.append(read)
                count_view[p_site - seg_start] += 1

        if do_warn == 1:
            warnings.warn("Data contains read alignments shorter (%s nt) than offset (%s nt). Ignoring." % (read_length,self.offset),
                          DataWarning)

        return reads_out, count_array


cdef class cVariableFivePrimeMapFactory(object):
    """Fiveprime-variable mapping for :py:meth:`BAMGenomeArray.set_mapping`.
    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
    The offset for a read of a given length is supplied in `offset_dict[readlen]`
     
    Parameters
    ----------
    offset_dict : dict
        Dictionary mapping read lengths to offsets that should be applied
        to reads of that length during mapping. A special key, `'default'` may
        be supplied to provide a default value for lengths not specifically
        enumerated in `offset_dict`
    """
    cdef int [10000] forward_offsets
    cdef int [10000] reverse_offsets

    def __cinit__(self, dict offset_dict not None):
        cdef int [:] fw_view = self.forward_offsets
        cdef int [:] rc_view = self.reverse_offsets
        cdef int offset
        fw_view [:] = -1
        rc_view [:] = 1

        # store offsets in c arrays for fast lookup
        # value for `default` in 0th place
        for read_length, offset in offset_dict.items():
            if read_length == "default":
                fw_view[0] = offset
                rc_view[0] = - offset - 1
            else:
                if offset >= read_length:
                    warnings.warn("Offset %s longer than read length %s. Ignoring %s-mers." % (offset, read_length, read_length),
                    DataWarning)
                    continue
                fw_view[read_length] = offset
                rc_view[read_length] = - offset - 1
         
    def __call__(self, list reads not None, c_roitools.GenomicSegment seg not None):
        """Returns reads covering a region, and a count vector mapping reads
        to specific positions in the region, mapping reads at possibly varying
        offsets from the 5' end of each read as determined by `self.offset_dict`
 
        Parameters
        ----------
        seg : |GenomicSegment|
            Region of interest
        
        reads : list of :py:class:`pysam.AlignedSegment`
            Reads to map
            
        Returns
        -------
        :py:class:`numpy.ndarray`
            Vector of counts at each position in `seg`
        """
        cdef:
            long seg_start = seg.start
            long seg_end = seg.end
            long seg_len = seg_end - seg_start
            list reads_out = []
            int warn_length = 0
            int [10000] offsets = self.forward_offsets
            int read_length, offset, p_site
            int empty_val = -1
            int do_warning = 0
            AlignedSegment read
            np.ndarray[LONG_t,ndim=1] count_array = np.zeros(seg_len,dtype=LONG)
            long [:] count_view = count_array
            long psite

        if seg.strand == reverse_strand:
            offsets = self.reverse_offsets
            empty_val = 1

        for read in reads:
            read_positions = <list>read.positions
            read_length = len(read_positions)
            offset = offsets[read_length]

            # check value is defined, if not use default
            if offset == empty_val:
                offset = offsets[0]
                # if defaut not defined, ignore read and warn
                if offset == empty_val:
                    do_warning = 1
                    continue

            p_site = read_positions[offset]
            if p_site >= seg_start and p_site < seg_end:
                reads_out.append(read)
                count_view[p_site - seg.start] += 1

        if do_warning == 1:
            warnings.warn("No offset for reads of length %s in offset dict. Ignoring %s-mers." % (warn_length, warn_length),
                          DataWarning)

        return reads_out, count_array


cdef class SizeFilterFactory(object):
    """Create a read-length filter can be applied at runtime to a |BAMGenomeArray|
    using ::meth:`BAMGenomeArray.add_filter`
    
    Parameters
    ----------
    min : float, optional
        Minimum read length to pass filter, inclusive (Default: `1`)
    
    max : float or numpy.inf, optional
        Maximum read length to pass filter, inclusive (Default: infinity)
    
    Returns
    -------
    function
    """
    # logically these would be ints but numpy.inf is a float
    cdef:
        float min_
        float max_

    def __cinit__(self,float min = 1, float max = np.inf):
        """Create a read-length filter can be applied at runtime to a |BAMGenomeArray|
        using ::meth:`BAMGenomeArray.add_filter`
        
        Parameters
        ----------
        min : float >= 1, optional
            Minimum read length to pass filter, inclusive (Default: `1`)
        
        max : float or numpy.inf >= `min`, optional
            Maximum read length to pass filter, inclusive (Default: infinity)
        
        Returns
        -------
        function

         .. note::
            
            `min` and `max` are typed as floats only to allow `max`
            to be `numpy.inf` by default. Fractional values and
            numbers less than 1 are nonsensical, and will raise
            :class:`ValueError`
        """
        if max < min:
            raise ValueError("Alignment size filter: max read length must be >= min read length")

        if min < 1:
            raise ValueError("Alignment size filter: min read length must be >= 1. Got %s" % min)

        self.min_ = min
        self.max_ = max
        
    def __call__(self,AlignedSegment read not None):
        cdef unsigned int my_length = len(read.positions)
        return my_length >= self.min_ and my_length <= self.max_
