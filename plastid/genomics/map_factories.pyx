"""Cython implementations of :term:`mapping rules <mapping rule>`.

See also
--------
Module documentation for :mod:`~plastid.genomics.genome_array`
    For a detailed discussion of :term:`mapping rules <mapping rule>`.
"""

import warnings
import numpy as np
cimport numpy as np
cimport cython

from pysam.calignmentfile cimport AlignedSegment
from plastid.genomics.c_common cimport forward_strand, reverse_strand, unstranded
from plastid.genomics.roitools cimport GenomicSegment
from plastid.util.services.exceptions import DataWarning
from plastid.util.io.filters import CommentReader
from plastid.util.scriptlib.argparsers import _parse_variable_offset_file


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


cdef class CenterMapFactory:
    """Read mapping tool for :meth:`BAMGenomeArray.set_mapping`.
    A user-specified number of bases is removed from each side of each
    read alignment, and the `N` remaining bases are each apportioned `1/N`
    of the read count, so that the entire read is counted once.

    Parameters
    ----------
    nibble : int >= 0, optional
        Number of bases to remove from each side of read (Default: `0`)
    """

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

    @cython.boundscheck(False) # valid because indices are explicitly checked
    @cython.cdivision(True) # we can do this because we explicitly check map_length > 0
    def __call__(self, list reads not None, GenomicSegment seg not None):
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
        cdef unsigned int nibble = self.nibble

        # count array we will return
        cdef np.ndarray[DOUBLE_t,ndim=1] count_array = np.zeros(seg_len,dtype=DOUBLE)
        cdef double [:] count_view = count_array

        cdef list reads_out = []
        cdef list read_positions
 
        cdef:
            AlignedSegment read 
            int coord, map_length
            unsigned int read_length, i
            DOUBLE_t val
       
        for read in reads:
            read_positions = <list>read.positions
            read_length = len(read_positions)
            map_length = read_length - 2*nibble
            if map_length < 0:
                do_warn = 1
                continue
            elif map_length > 0:
                val = 1.0 / map_length
                for i in range(nibble,read_length - nibble):
                    coord = <long>(read_positions[i]) - seg_start
                    if coord >= 0 and coord < seg_len:
                        count_view[coord] += val

                reads_out.append(read)

        if do_warn == 1:
            warnings.warn("Data contains read alignments shorter than `2*nibble` value of %s nt. Ignoring these." % 2*nibble,
                  DataWarning)
        
        return reads_out, count_array

    property nibble:
        """Number of positions to trim from each side of read alignment before assigning genomic positions."""
        def __get__(self):
            return self.nibble
        def __set__(self, unsigned int val):
            self.nibble = val


cdef class FivePrimeMapFactory:
    """Fiveprime mapping factory for :py:meth:`BAMGenomeArray.set_mapping`.
    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
     
    Parameters
    ----------
    offset : int >= 0, optional
        Offset from 5' end of read, in direction of 3' end, at which 
        reads should be counted (Default: `0`)
    """

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

    def __call__(self, list reads not None, GenomicSegment seg not None):
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

        if seg.c_strand == reverse_strand:
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

    property offset:
        """Distance from 5' end of read at which to assign reads"""
        def __get__(self):
            return self.offset
        def __set__(self, int val):
            self.offset = val

     
cdef class ThreePrimeMapFactory:
    """Threeprime mapping factory for :py:meth:`BAMGenomeArray.set_mapping`.
    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
     
    Parameters
    ----------
    offset : int >= 0, optional
        Offset from 3' end of read, in direction of 5' end, at which 
        reads should be counted (Default: `0`)
    """

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

    def __call__(self, list reads not None, GenomicSegment seg not None):
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

        if seg.c_strand != reverse_strand:
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


    property offset:
        """Distance from 3' end of read at which to assign reads"""
        def __get__(self):
            return self.offset
        def __set__(self, int val):
            self.offset = val


cdef class VariableFivePrimeMapFactory:
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

    def __cinit__(self, dict offset_dict not None):
        cdef:
            int [:] fw_view = self.forward_offsets
            int [:] rc_view = self.reverse_offsets
            int offset
            object read_length # can be int or str
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
    
    @staticmethod
    def from_file(object fn_or_fh):
        """Create a :class:`VariableFivePrimeMapFactory` from a text file.

        The text file should be formatted as by the :mod:`~plastid.bin.psite`
        script: a tab-delimited two-column text file in which the first column
        is a read length (or the special word 'default'), and the second column
        an offset from the 5' end of the read corresponding to that read length.
        For example::

            25      12
            26      12
            27      13
            28      13
            29      14
            30      14
            default 14


        Parameters
        ----------
        fn_or_fh : str or open file-like
            If str, it should be a path to the text file. If file-like, it 
            should be an open file-like object pointing to the data that
            would be in the text file.
        """
        if isinstance(fn_or_fh,str):
            fh = CommentReader(open(str))
        else:
            fh = CommentReader(fn_or_fh)

        return VariableFivePrimeMapFactory(_parse_variable_offset_file(fh))
         
    def __call__(self, list reads not None, GenomicSegment seg not None):
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
            int no_offset_length = 0
            int bad_offset_length = 0
            int [10000] offsets = self.forward_offsets
            int empty_val = -1
            int do_no_offset_warning = 0
            int do_bad_default_warning = 0
            int read_length, offset, p_site
            AlignedSegment read
            np.ndarray[LONG_t,ndim=1] count_array = np.zeros(seg_len,dtype=LONG)
            long [:] count_view = count_array
            long psite

        if seg.c_strand == reverse_strand:
            offsets = self.reverse_offsets
            empty_val = 1

        for read in reads:
            read_positions = <list>read.positions
            read_length = len(read_positions)
            fo = self.forward_offsets[read_length] # needed for read length test
            offset = offsets[read_length]

            # check value is defined, if not use default
            if offset == empty_val:
                offset = offsets[0]
                fo = self.forward_offsets[0]
                # if defaut not defined, ignore read and warn
                if offset == empty_val:
                    do_no_offset_warning = 1
                    no_offset_length = read_length
                    continue
            
            if fo >= read_length:
                do_bad_default_warning = 1
                bad_offset_length = read_length
                continue

            p_site = read_positions[offset]
            if p_site >= seg_start and p_site < seg_end:
                reads_out.append(read)
                count_view[p_site - seg.start] += 1

        if do_no_offset_warning == 1:
            warnings.warn("No offset for reads of length %s nt in offset dict. Ignoring these." % (no_offset_length),
                          DataWarning)

        if do_bad_default_warning == 1:
            warnings.warn("'Default' offset (%s nt) exceeds length of some reads found (%s nt). Ignoring these." % (self.forward_offsets[0], bad_offset_length),
                          DataWarning)

        return reads_out, count_array


cdef class SizeFilterFactory:
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

    def __cinit__(self,int min = 1, int max = -1):
        """Create a read-length filter can be applied at runtime to a |BAMGenomeArray|
        using ::meth:`BAMGenomeArray.add_filter`
        
        Parameters
        ----------
        min : int >= 1, optional
            Minimum read length to pass filter, inclusive (Default: `1`)
        
        max : int, optional
            Maximum read length to pass filter, inclusive. If `-1`, then
            there is no maximum length filter. (Default: `-1`, no filter.)
        
        Returns
        -------
        function
        """
        if max != -1 and max < min:
            raise ValueError("Alignment size filter: max read length must be >= min read length")

        if min < 1:
            raise ValueError("Alignment size filter: min read length must be >= 1. Got %s" % min)

        self.min_ = min
        self.max_ = max
        
    def __call__(self,AlignedSegment read not None):
        cdef unsigned int my_length = len(read.positions)
        return my_length >= self.min_ and (my_length <= self.max_ or self.max_ == -1)

