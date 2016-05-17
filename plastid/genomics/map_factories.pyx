""":term:`Mapping functions` used by |BAMGenomeArray|

.. contents::
   :local:

Summary
-------
Numerous sequencing assays encode interesting biology in various properties of
read alignments. For example:

  - Ribosome profiling encodes the positions of ribosomal P-sites as a
    function of both read length and alignment coordinates
 
  - Bisulfite sequencing encodes methylation status as C-to-T transitions
    within read alignments

  - DMS-Seq and Pseudouridine profiling encode unstructured regions of RNA or
    sites of pseudouridine modification, respectively, in the 5' termini fo
    their read aligments

Extended summary
----------------

plastid uses configurable mapping functions to
decode the biology of interest from read alignments.

This design enables plastid's tools to operate on sequencing data from
virtually any NGS assay, provided the appropriate  mapping function and
parameters. Mapping rules are described in the following section: 

This module contains Cython implementations of mapping functions for BAM
files. At present, the following are included:


*Fiveprime end mapping*
   Each read alignment is mapped to its 5' end, or at a fixed offset (in
   nucleotides) from its 5' end
        
*Variable fiveprime end mapping*
   Each read alignment is mapped at a fixed distance from its 5' end, where
   the distance is determined by the length of the read alignment.
     
   This is used for :term:`ribosome profiling` of yeast (:cite:`Ingolia2009`)
   and mammalian cells (:cite:`Ingolia2011`).

*Stratified variable fiveprime end mapping*
   This multidimensional mapping rule behaves as variable fiveprime end mapping,
   and additionally returns arrays that are stratified by read length.
   
   For each region of interest, a 2D array is returned, in which each row
   represents a given read length, and each column a nucleotide position in the
   region of interest. In other words, summing over the columns produces the
   same array that would be given by variable fiveprime end mapping.

*Threeprime end mapping*
   Each read alignment is mapped to its 3' end, or at a fixed
   offset (in nucleotides) from its 3' end.
    
*Entire* or *Center-weighted mapping*
   Zero or more positions are trimmed from each end of the read alignment,
   and the remaining `N` positions in the alignment are incremented by `1/N`
   read counts (so that each read is still counted once, when integrated
   over its mapped length).
     
   This is also used for ribosome profiling of E. coli and D. melanogaster, and
   RNA-seq. 



Examples
--------

Mapping functions are rarely called directly. Instead, they are invoked via
:meth:`plastid.genomics.genome_array.BAMGenomeArray.set_mapping`:

 .. code-block:: python

    # open BAM file
    >>> ga = BAMGenomeArray(["some_file.bam"])

    # set mapping 15 bases from 3' end of read
    >>> ga.set_mapping(ThreePrimeMapFactory(offset=15))

    # set mapping to 5' end of read
    >>> ga.set_mapping(FivePrimeMapFactory())

    # cut 3 nt off each end of the read, and map fractionally over the rest
    >>> ga.set_mapping(CenterMapFactory(nibble=12))

    # access reads in some region
    >>> ga[GenomicSegment("chrI",10000,100000,"+")]
    # output is numpy array

Other mapping functions are similarly invoked. See their docstrings for further details.


Implementation
--------------

Mapping functions must be callables, and can therefore be implemented as functions
or classes. In order to give them configurable parameters, one could create
function factories or use :func:`functools.partial` to supply the extra
parameters. Here, we created classes that produce callable instances that mimic
functions. We did this because early on we ran into trouble with limitations of
:mod:`pickle` and :mod:`multiprocessing`, which can only capture functions and
classes defined in the top-level scope. This limitation is no longer relevant
to our purposes but is preserved here in history.

Generally, the callable must accept two arguments: a list of read alignments
(as :class:`pysam.AlignedSegment`) and a |GenomicSegment| representing a query
interval. The callable must then return two values: a list of read alignments,
and a :class:`numpy.ndarray` of values.

Generally, all reads are mapped within the callable. Only those that contribute
data to coordinates contained within the query interval are returned.

The :class:`numpy.ndarray` of values may be multidimensional. The only stipulation
is that the last dimension represent the positions in the query interval. For
a 2D array, the positions would thus be columns. See
|StratifiedVariableFivePrimeMapFactory| below for an example.



See also
--------
Module documentation for :mod:`~plastid.genomics.genome_array`
    For a detailed discussion of :term:`mapping functions <mapping rule>`.
"""

import numpy as np
cimport numpy as np
cimport cython

from pysam.calignmentfile cimport AlignedSegment
from plastid.genomics.c_common cimport forward_strand, reverse_strand, unstranded
from plastid.genomics.roitools cimport GenomicSegment
from plastid.util.services.exceptions import DataWarning, warn, warn_onceperfamily
from plastid.util.io.filters import CommentReader
from plastid.util.scriptlib.argparsers import _parse_variable_offset_file

DEF _BAD_OFFSET = -1
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
    """
    CenterMapFactory(nibble=0)

    Read mapping tool for :meth:`BAMGenomeArray.set_mapping`.
    `nibble` positions  is removed from each side of each read alignment, and
    the `N` remaining positions are each apportioned `1/N` of the read count,
    so that the entire read is counted once.

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
        list
            List of reads that were mapped into the output array

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
            warn_onceperfamily("Data contains read alignments shorter than `2*nibble` value of '%s' nt. Ignoring these." % (2*nibble),
                  DataWarning)
        
        return reads_out, count_array

    property nibble:
        """Number of positions to trim from each side of read alignment before assigning genomic positions."""
        def __get__(self):
            return self.nibble
        def __set__(self, unsigned int val):
            self.nibble = val


cdef class FivePrimeMapFactory:
    """
    FivePrimeMapFactory(offset=0)
    
    Fiveprime mapping factory for :py:meth:`BAMGenomeArray.set_mapping`.
    Reads are mapped at `offset` nucleotides from the fiveprime end of their
    alignment.
     
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
            warn_onceperfamily("Data contains read alignments shorter than offset (%s nt). Ignoring." % (self.offset),
                 DataWarning)

        return reads_out, count_array

    property offset:
        """Distance from 5' end of read at which to assign reads"""
        def __get__(self):
            return self.offset
        def __set__(self, int val):
            self.offset = val

     
cdef class ThreePrimeMapFactory:
    """
    ThreePrimeMapFactory(offset=0)
    
    Threeprime mapping factory for :py:meth:`BAMGenomeArray.set_mapping`.
    Reads are mapped `offset` nucleotides from the threeprime end of their alignments.
     
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
            warn_onceperfamily("Data contains read alignments shorter than offset (%s nt). Ignoring." % self.offset,
                 DataWarning)

        return reads_out, count_array


    property offset:
        """Distance from 3' end of read at which to assign reads"""
        def __get__(self):
            return self.offset
        def __set__(self, int val):
            self.offset = val


cdef class VariableFivePrimeMapFactory:
    """
    VariableFivePrimeMapFactory(offset_dict)
    
    Fiveprime-variable mapping for :py:meth:`BAMGenomeArray.set_mapping`.
    Reads are mapped at a user-specified offsets from the fiveprime end of each
    alignment. The offset for a read of a given length is supplied in `offset_dict[readlen]`
     
    Parameters
    ----------
    offset_dict : dict
        Dictionary mapping read lengths to offsets that should be applied
        to reads of that length during mapping. A special key, `'default'` may
        be supplied to provide a default value for lengths not specifically
        enumerated in `offset_dict`
    """

    def __cinit__(self, dict offset_dict, *args):
        """
        Parameters
        ----------
        offset_dict : dict
            Dictionary mapping read lengths to offsets that should be applied
            to reads of that length during mapping. A special key, `'default'` may
            be supplied to provide a default value for lengths not specifically
            enumerated in `offset_dict`
        """
        cdef:
            int [:] fw_view = self.forward_offsets
            int [:] rc_view = self.reverse_offsets
            int offset
            object read_length # can be int or str
            int i
            
        fw_view [:] = _BAD_OFFSET
        rc_view [:] = _BAD_OFFSET
        
        # this might seem silly, because nobody would every make a 
        # FivePrimeVariableMapFactory without an offset_dict, but this allows
        # reuse of this __cinit__ in StratifiedVariableFivePrimeMapFactory.
        if offset_dict is None:
            offset_dict = { "default" : 0 }
            
        if "default" in offset_dict:
            default = int(offset_dict["default"])
            fw_view[default + 1:] = default
            i = default + 1
            while i < len(rc_view):
                rc_view[i] = i - default - 1
                i += 1

        # store offsets in c arrays for fast lookup
        # value for `default` in 0th place
        for read_length, offset in offset_dict.items():
            if read_length != "default":
                if offset >= read_length:
                    if read_length >= default:
                        warn_onceperfamily("Given offset '%s' longer than read length '%s'. Falling back to default '%s'." % (offset, read_length, default),
                             DataWarning)
                        offset = default
                    else:
                        warn_onceperfamily("Given offset '%s' and default '%s' are longer than read length '%s'. Ignoring %s-mers." % (offset, default, read_length, read_length),
                             DataWarning)
                    continue

                fw_view[read_length] = offset
                rc_view[read_length] = read_length - offset - 1
    
    @staticmethod
    def from_file(object fn_or_fh):
        """
        from_file(object fn_or_fh)
        
        Create a :class:`VariableFivePrimeMapFactory` from a text file.

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

    @cython.boundscheck(False) # valid because indices are explicitly checked in __cinit__
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
        list
            List of reads that were mapped into the output array

        :py:class:`numpy.ndarray`
            Vector of counts at each position in `seg`
        """
        cdef:
            int [:] offsets = self.forward_offsets
            
            long seg_start = seg.start
            long seg_end   = seg.end
            long seg_len   = seg_end - seg_start
            
            int no_offset_length
            int do_no_offset_warning = 0
            
            int            read_length, offset
            long           p_site
            AlignedSegment read
            
            np.ndarray count_array = np.zeros(seg_len,dtype=LONG)
            long [:]   count_view  = count_array
            list read_positions
            list reads_out       = []

        if seg.c_strand == reverse_strand:
            offsets = self.reverse_offsets

        for read in reads:
            read_positions = <list>read.positions
            read_length    = len(read_positions)
            offset         = offsets[read_length]

            if offset == _BAD_OFFSET:
                do_no_offset_warning = 1
                no_offset_length = read_length
                continue

            p_site = read_positions[offset]
            if p_site >= seg_start and p_site < seg_end:
                reads_out.append(read)
                count_view[p_site - seg.start] += 1

        if do_no_offset_warning == 1:
            warn_onceperfamily("No usable offset for reads of length %s nt in offset dict. Ignoring these." % (no_offset_length),
                 DataWarning)

        return reads_out, count_array


cdef class StratifiedVariableFivePrimeMapFactory(VariableFivePrimeMapFactory):
    """
    StratifiedVariableFivePrimeMapFactory(offset_dict, min = 25, max = 35)
    
    Fiveprime-variable mapping for :py:meth:`BAMGenomeArray.set_mapping`, stratified by read length.
    
    Reads are mapped into a 2D array of counts at each position in a region of
    interest (column), stratified by read length (row). Mapping is performed
    by applying a user-specified offset from the fiveprime end of each alignment.
    A unique offset may be supplied for each read length. 
    
    Reads outside of pre-specified minimum and maximum lengths are ignored.


    Parameters
    ----------
    offset_dict : dict or `None`
        Dictionary mapping read lengths to offsets that should be applied
        to reads of that length during mapping. A special key, `'default'` may
        be supplied to provide a default value for lengths not specifically
        enumerated in `offset_dict.` If `None`, all reads are mapped to their
        5' ends.
        
    min_length : int
        Minimum length read to map
        
    max_length : int
        Maximum length read to map, inclusive

    
    Attributes
    ----------
    row_keys : numpy.ndarray
        1D numpy array indicating which read length is assigned to which row

    shape : list of integers
        size of each axis of output, excluding the final axis, which is
        specified the length of incoming GenomicSegments.
        
     .. warning::
     
        The :meth:`~plastid.genomics.genome_array.BAMGenomeArray.to_variable_step`,
        :meth:`~plastid.genomics.genome_array.BAMGenomeArray.to_bedgraph`,
        :meth:`~plastid.genomics.genome_array.BAMGenomeArray.to_genome_array`,
        of :class:`BAMGenoemArray` do not support multi-dimensional map
        factories.-
    """
    def __cinit__(self, dict offset_dict, int min = 25, int max = 35):
        """
        Parameters
        ----------
        offset_dict : dict
            Dictionary mapping read lengths to offsets that should be applied
            to reads of that length during mapping. A special key, `'default'` may
            be supplied to provide a default value for lengths not specifically
            enumerated in `offset_dict`
            
        min_length : int
            Minimum length read to map
            
        max_length : int
            Maximum length read to map, inclusive
        """
        if max <= min:
            raise ValueError("Max length '%s' must be >= min length '%s'. " % (max,min))
        
        self.min_length  = min
        self.max_length  = max
        self._numlengths = max - min + 1

    @cython.boundscheck(False) # valid because indices are explicitly checked in VariableFivePrimeMapFactory.__cinit__
    def __call__(self, list reads not None, GenomicSegment seg not None):
        """
        Map reads covering `seg` into a 2D count array, in which reads at
        each nucleotide position (column) are stratified by read length (row),
        and, optionally, offset from their 5' ends.
 
        Parameters
        ----------
        seg : |GenomicSegment|
            Region of interest
        
        reads : list of :py:class:`pysam.AlignedSegment`
            Reads to map
            
        Returns
        -------
        list
            List of reads that were mapped into the output array

        :py:class:`numpy.ndarray`
            2D array, in which each cell value is the number of counts
            of a given length (row) at  given position (column)
        """
        cdef:
            long seg_start             = seg.start
            long seg_end               = seg.end
            long seg_len               = seg_end - seg_start

            int [10000] offsets        = self.forward_offsets
            int num_lengths            = self._numlengths
            int min_length             = self.min_length
            
            list reads_out             = []
            int read_length, offset, p_site
            
            AlignedSegment read
            
            np.ndarray count_array = np.zeros((num_lengths,seg_len),dtype=LONG)
            long [:,:] count_view  = count_array
            long psite

        if seg.c_strand == reverse_strand:
            offsets = self.reverse_offsets

        for read in reads:
            read_positions = <list>read.positions
            read_length = len(read_positions)
            if read_length >= min_length and read_length <= self.max_length:

                offset = offsets[read_length]
                p_site = read_positions[offset]
                
                if p_site >= seg_start and p_site < seg_end:
                    reads_out.append(read)
                    count_view[read_length - min_length, p_site - seg_start] += 1

        return reads_out, count_array

    property row_keys:
        """numpy array of read lengths corresponding to each row of mapped data."""
        def __get__(self):
            return np.arange(self.min_length,self.max_length+1)
        
    property shape:
        """list containing size of each axis of output, excluding the final axis,
        which is specified the length of incoming GenomicSegments."""
        def __get__(self):
            return [self._numlengths]


cdef class SizeFilterFactory:
    """
    SizeFilterFactory(min = 1, max = -1)
    
    Create a read-length filter can be applied at runtime to a |BAMGenomeArray|
    using ::meth:`BAMGenomeArray.add_filter`
    
    Parameters
    ----------
    min : float, optional
        Minimum read length to pass filter, inclusive (Default: `1`)
    
    max : float or numpy.inf, optional
        Maximum read length to pass filter, inclusive. If set to `-1`,
        then there is no maximum length filter. (Default: -1, no filter)
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

