#!/usr/bin/env python
"""GenomeArrays are randomly-accessible array-like structures that map quantitative data
or :term:`read alignments` to genomic positions. In this way, they function like
`numpy`_ arrays, except over coordinates specified by |GenomicSegments|
rather than slices of `start:end`. For example, to fetch data covering nucleotdies
2000-2500 on the plus strand of chromosome `'chrA'`::

    >>> genome_array = GenomeArray()
    >>> genome_array.add_from_wiggle(open("some_file_fw.wig"),"+")
    >>> genome_array.add_from_wiggle(open("some_file_rc.wig"),"-") 
    >>> segment = GenomicSegment('chrA',2000,2500,'+')
    >>> my_data = some_genome_array[segment]
        # some numpy array


More commonly, GenomeArrays are used with |SegmentChains| and |Transcripts|, which
via their :meth:`~plastid.genomics.roitools.SegmentChain.get_counts` methods
can fetch arrays of data in which each array position corresponds to a position,
moving from 5' to 3', in the spliced |SegmentChain| or |Transcript|:: 

    >>> seg1 = GenomicSegment('chrA',2000,2500,'-')
    >>> seg2 = GenomicSegment('chrA',3000,3250','-')
    >>> seg3 = GenomicSegment('chrA',4000,4200,'-')
    >>> segment_chain = SegmentChain(seg1,seg2,seg3,ID="my_chain")
    >>> my_data = segment_chain.get_counts(some_genome_array)
        # some numpy array
        

Types of GenomeArrays
---------------------

Several subclasses are provided, depending on how the :term:`read alignments`
or quantitative data is stored. All GenomeArrays present the same interfaces for:

    - retrieving vectors (as :py:class:`numpy.ndarray` s) of data at each position
      of a |GenomicSegment|
      
    - tabulating the sum of the :term:`read alignments` or :term:`counts`
      in a dataset via their :meth:`~AbstractGenomeArray.sum` methods
    
    - toggling normalization of fetched :term:`counts` to reads-per-million
      via their :meth:`~AbstractGenomeArray.set_normalize` methods

    - export of data to variableStep `Wiggle`_ and `bedGraph`_ formats via their
      methods :meth:`~GenomeArray.to_variable_step` and
      :meth:`~GenomeArray.to_bedgraph`, respectively
    

The following subclasses of GenomeArrays provide their own features in addition 
to those above:

    |GenomeArray|
        The most basic GenomeArray. It allows:
        
            - setting values within the array
    
            - import of quantitative data from `Wiggle`_ and `bedGraph`_ files
            
            - import of :term:`read alignments` from native `bowtie`_ alignments,
              and their convertion to :term:`counts` under configurable
              :term:`mapping rules <mapping rule>`.
    
            - mathematical operations such as addition and subtraction of scalars,
              or positionwise addition and subtraction of other |GenomeArray| s
    
    
    |SparseGenomeArray|
        A slower but more memory-efficient implementation of |GenomeArray|,
        useful for large genomes or on computers with limited memory.


    |BAMGenomeArray|
        A GenomeArray for `BAM` files. Fetches alignments from one or more
        `BAM`_ files, and applies :term:`mapping rules <mapping rule>` to
        convert these alignments to to vectors of :term:`counts`.
        
        `BAM`_ files give |BAMGenomeArray| several advantages over |GenomeArray|:
        
            - |BAMGenomeArray| objects can be much faster to create and far more
              memory-efficient than |GenomeArray| or |SparseGenomeArray|.
            
            - :term:`Mapping rules <mapping rule>` and filter functions may be
              changed at runtime, with no speed cost, rather than having to be
              set at import time (as is the case for import of `bowtie`_ files
              into a |GenomeArray|). See :meth:`BAMGenomeArray.set_mapping` for details.
                
            - All of the rich properties of :term:`read alignments` that are
              preserved in `BAM`_ files (e.g. alignment mismatches, read lengths,
              et c) are accessible, and may be used by :term:`mapping rules <mapping rule>`
              or filter functions. 
              
            - If necessary, a |BAMGenomeArray| can be converted to a |GenomeArray|
        
        
        However, because a |BAMGenomeArray| is a view of the data in an
        underlying `BAM`_ file, |BAMGenomeArrays| do not support
        *setter* operations. Instead, mathematical operations, if desired, may be:
        
            - applied to vectors of counts, once fetched from the |BAMGenomeArray|
            
            - performed on a |GenomeArray| or |SparseGenomeArray| created from
              a |BAMGenomeArray| via :meth:`BAMGenomeArray.to_genome_array`


Mapping rules
-------------
In order to convert :term:`read alignments` to :term:`counts`, a :term:`mapping rule`
must be applied. Depending on the data type and the purpose of the experiment,
different rules may be appropriate. The following mapping rules are provided,
although users are encouraged to provide their own if needed. These include:

    #.  *Fiveprime end mapping:*
        Each read alignment is mapped to its 5' end, or at a fixed
        distance from its 5' end. This is common for RNA-seq or CLiP-seq
        experiments.
        
    #.  *Variable fiveprime end mapping:*
        Each read alignment is mapped at a fixed distance from its
        5' end, where the distance is determined by the length of the read
        alignment. This is commonly used for mapping :term:`ribosome-protected footprints`
        to their P-sites in :term:`ribosome profiling` experiments.
    
    #.  *Threeprime end mapping:*
        Each read alignment is mapped to its 3' end, or at a fixed
        distance from its 3' end.
    
    #.  *Entire,* *Center-weighted,* or *nibble mapping:*
        Zero or more positions are trimmed from each end of the read alignment,
        and the remaining `N` positions in the alignment are incremented by `1/N`
        read counts (so that each read is still counted once, when integrated
        over its mapped length). This is also occasionally used for 
        :term:`ribosome profiling` or RNA-seq. 

Implementations of each mapping rule are provided. For |GenomeArray| and
|SparseGenomeArray|, functions are provided as arguments to the
method :meth:`~GenomeArray.add_from_bowtie`. For |BAMGenomeArray|, a function
is generated from one of the :term:`factories <factory>` below, and passed to
the method :meth:`BAMGenomeArray.set_mapping`:

======================   ====================================    =======================================
**Mapping rule**         |GenomeArray|, |SparseGenomeArray|      |BAMGenomeArray|
----------------------   ------------------------------------    ---------------------------------------

Fiveprime                :func:`five_prime_map`                  :py:func:`FivePrimeMapFactory`

Fiveprime variable       :func:`variable_five_prime_map`         :py:func:`VariableFivePrimeMapFactory`

Threeprime               :func:`three_prime_map`                 :py:func:`ThreePrimeMapFactory`

Center/entire/nibble     :func:`center_map`                      :py:func:`CenterMapFactory`
======================   ====================================    =======================================


"""
__date__ =  "May 3, 2011"
__author__ = "joshua"
from abc import abstractmethod
import itertools
import operator
import copy
import warnings
import numpy
import scipy.sparse
from collections import OrderedDict
from plastid.readers.wiggle import WiggleReader
from plastid.readers.bowtie import BowtieReader
from plastid.genomics.roitools import GenomicSegment, SegmentChain
from plastid.util.services.mini2to3 import xrange, ifilter
from plastid.util.services.exceptions import DataWarning

from plastid.genomics.map_factories import *

MIN_CHR_SIZE = 10*1e6 # 10 Mb minimum size for unspecified chromosomes 

#
##===============================================================================
## Factories for mapping functions for BAMGenomeArray or other structures
## Each factory returns a function that takes a list of pysam.AlignedSegments
## and a GenomicSegment, and returns a list of pysam.AlignedSegments that mapped
## to that interval under the mapping rules specified in that function, and
## a corresponding vector of read counts mapping to each nucleotide position
## in the GenomicSegment
##===============================================================================
#
#def CenterMapFactory(nibble=0):
#    """Returns a center-mapping function for :py:meth:`BAMGenomeArray.set_mapping`.
#    A user-specified number of bases is removed from each side of each
#    read alignment, and the `N` remaining bases are each apportioned `1/N`
#    of the read count, so that the entire read is counted once.
#     
#    Parameters
#    ----------
#    nibble : int, optional
#        Number of bases to remove from each side of read (Default: `0`)
#    
#    Returns
#    -------
#    function
#        Mapping function
#    """
#    # docstring of function we will return.
#    # set it up here so we can put the string substitution in it
#    docstring = """Returns reads covering a region, and a count vector mapping reads
#        to specific positions in the region. %s bases are trimmed from each
#        side of the read, and each of the `N` remaining alignment positions
#        are incremented by `1/N`
#        
#        Parameters
#        ----------
#        seg : |GenomicSegment|
#            Region of interest
#        
#        reads : list of :py:class:`pysam.AlignedSegment`
#            Reads to map
#            
#        Returns
#        -------
#        :py:class:`numpy.ndarray`
#            Vector of counts at each position in `seg`
#        """ % nibble
# 
#    def map_func(reads,seg):
#        seg_start = seg.start
#        seg_end = seg.end
#        seg_positions = frozenset(range(seg_start,seg_end))
#        reads_out = []
#        count_array = numpy.zeros(seg_end - seg_start,dtype=float)
#
#        for read in reads:
#            read_positions = read.positions
#            read_length = len(read_positions)
#            if read_length <= 2*nibble:
#                warnings.warn("File contains read alignments shorter (%s nt) than `2*'nibble'` value of %s nt. Ignoring these." % (read_length,2*nibble),
#                              DataWarning)
#                continue
#
#            if nibble > 0:
#                read_positions = read_positions[nibble:-nibble]
#
#            overlap = numpy.array(list(frozenset(read_positions) & seg_positions),dtype=int)
#            if len(overlap) > 0:
#                overlap -= seg_start
#                reads_out.append(read)
#                val = 1.0 / (read_length - 2*nibble)
#                count_array[overlap] += val
#
#        return reads_out, count_array
#     
#    map_func.__doc__ = docstring
#    map_func.__mapping__ = "center"
#    return map_func
# 
# 
#def FivePrimeMapFactory(offset=0):
#    """Returns a fiveprime mapping function for :py:meth:`BAMGenomeArray.set_mapping`.
#    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
#     
#    Parameters
#    ----------
#    offset : int, optional
#        Offset from 5' end of read, in direction of 3' end, at which 
#        reads should be counted (Default: `0`)
#    
#    Returns
#    -------
#    function
#        Mapping function
#    """
#     
#    # docstring of function we will return.
#    # set it up here so we can put the string substitution in it
#    docstring = """Returns reads covering a region, and a count vector mapping reads
#        to specific positions in the region, mapping reads at %s
#        from the fiveprime end of each read.
# 
#        Parameters
#        ----------
#        reads : list of :py:class:`pysam.AlignedSegment`
#            Reads to map
#            
#        seg : |GenomicSegment|
#            Region of interest
#
#            
#        Returns
#        -------
#        list
#            List of :py:class:`pysam.AlignedSegment` that map to region
#            
#        :py:class:`numpy.ndarray`
#            Vector of counts at each position in `seg`
#        """ % offset
#         
#    def map_func(reads,seg):
#        reads_out = []         
#        seg_start = seg.start
#        seg_end = seg.end
#        seg_len = seg_end - seg_start
#        count_array = numpy.zeros(seg_len)
#
#        read_offset = offset
#        if seg.strand == "-":
#            read_offset = - offset - 1
#
#        for read in reads:
#            read_positions = read.positions
#            read_length = len(read_positions)
#            if abs(offset) >= read_length:
#                warnings.warn("File contains read alignments shorter (%s nt) than offset (%s nt). Ignoring." % (read_length,offset),
#                              DataWarning)
#                continue
#            
#            p_site = read_positions[read_offset]
#            if p_site >= seg_start and p_site < seg_end:
#                reads_out.append(read)
#                count_array[p_site - seg_start] += 1
#        return reads_out, count_array
#     
#    map_func.__doc__ = docstring
#    map_func.__mapping__ = "fiveprime"    
#    return map_func
# 
# 
#def ThreePrimeMapFactory(offset=0):
#    """Returns a threeprime mapping function for :py:meth:`BAMGenomeArray.set_mapping`.
#    Reads are mapped at a user-specified offset from the 3' end of the alignment,
#    in the direction of the 5' end
#     
#    Parameters
#    ----------
#    offset : int, optional
#        Offset from 3' end of read, in direction of 5' end, at which 
#        reads should be counted (Default: *0*)
#    
#    Returns
#    -------
#    function
#        Mapping function
#    """
#     
#    # docstring of function we will return.
#    # set it up here so we can put the string substitution in it
#    docstring = """Returns reads covering a region, and a count vector mapping reads
#        to specific positions in the region, mapping reads at %s
#        from the threeprime end of each read.
# 
#        Parameters
#        ----------
#        seg : |GenomicSegment|
#            Region of interest
#        
#        reads : list of :py:class:`pysam.AlignedSegment`
#            Reads to map
#            
#        Returns
#        -------
#        :py:class:`numpy.ndarray`
#            Vector of counts at each position in `seg`
#        """ % offset
#    def map_func(reads,seg):
#        seg_start = seg.start
#        seg_end   = seg.end
#        seg_len   = seg_end - seg_start
#        reads_out = []
#        count_array = numpy.zeros(seg_len)
#
#        read_offset = offset
#        if seg.strand == "+":
#            read_offset = - offset - 1
#
#        for read in reads:
#            read_positions = read.positions
#            read_length = len(read_positions)
#            if offset >= read_length: 
#                warnings.warn("File contains read alignments shorter (%s nt) than offset (%s nt). Ignoring." % (read_length,offset),
#                              DataWarning)
#                continue
#
#            p_site = read.positions[read_offset]  
#            if p_site >= seg_start and p_site < seg_end:
#                reads_out.append(read)
#                count_array[p_site - seg_start] += 1
#
#        return reads_out, count_array
#     
#    map_func.__doc__ = docstring
#    map_func.__mapping__ = "threeprime"    
#    return map_func
# 
## JGD modified 2013-12-16
## APF 2013-11-18
#def VariableFivePrimeMapFactory(offset_dict):
#    """Returns a fiveprime-variable mapping function for :py:meth:`BAMGenomeArray.set_mapping`.
#    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
#    The offset for a read of a given length is supplied in `offset_dict[readlen]`
#     
#    Parameters
#    ----------
#    offset_dict : dict
#        Dictionary mapping read lengths to offsets that should be applied
#        to reads of that length during mapping. A special key, `'default'` may
#        be supplied to provide a default value for lengths not specifically
#        enumerated in `offset_dict`
#    
#    Returns
#    -------
#    function
#        Mapping function    
#    """
#     
#    # docstring of function we will return.
#    docstring = """Returns reads covering a region, and a count vector mapping reads
#        to specific positions in the region, mapping reads at possibly varying
#        offsets from the 5' end of each read.
# 
#        Parameters
#        ----------
#        seg : |GenomicSegment|
#            Region of interest
#        
#        reads : list of :py:class:`pysam.AlignedSegment`
#            Reads to map
#            
#        Returns
#        -------
#        :py:class:`numpy.ndarray`
#            Vector of counts at each position in `seg`
#        """
#         
#    def map_func(reads,seg,offset_dict=offset_dict):
#        reads_out = []         
#        count_array = numpy.zeros(len(seg))
#        for read in reads:
#            read_length = len(read.positions)
#            if read_length in offset_dict:
#                offset = offset_dict[read_length]
#            elif "default" in offset_dict:
#                offset = offset_dict["default"]
#            else:
#                warnings.warn("No offset for reads of length %s in offset dict. Ignoring %s-mers." % (read_length, read_length),
#                              DataWarning)
#                continue
#
#            if offset >= read_length:
#                warnings.warn("Offset %s longer than read length %s. Ignoring %s-mers." % (offset, read_length, read_length),
#                              DataWarning)
#                continue
#
#            if seg.strand == "+":
#                p_site = read.positions[offset]
#            else:
#                p_site = read.positions[-offset - 1]
#             
#            if p_site >= seg.start and p_site < seg.end:
#                reads_out.append(read)
#                count_array[p_site - seg.start] += 1
#
#        return reads_out, count_array
#     
#    map_func.__doc__ = docstring
#    map_func.__mapping__ = "fiveprime_variable"    
#    return map_func
# 
# 
#
##===============================================================================
## Factory functions to filter reads for BAMGenomeArrays
##===============================================================================
#
#def SizeFilterFactory(min=1,max=numpy.inf):
#    """Return a function to filter reads based on length, which can be applied
#    at runtime to a |BAMGenomeArray| using :py:meth:`BAMGenomeArray.add_filter`
#    
#    Parameters
#    ----------
#    min : int, optional
#        Minimum read length to pass filter, inclusive (Default: `1`)
#    
#    max : int or numpy.inf, optional
#        Maximum read length to pass filter, inclusive (Default: infinity)
#    
#    Returns
#    -------
#    function
#    """
#    if max < min:
#        raise ValueError("Alignment size filter: max read length must be >= min read length")
#    
#    def my_func(x):
#        my_length = len(x.positions)
#        return my_length >= min and my_length <= max
#    
#    return my_func
#
#

#===============================================================================
# INDEX: Mapping functions for GenomeArray and SparseGenomeArray.
#        these are used by add_from_bowtie() and add_from_tagalign()
#        to map read alignments to specific sites.
#
#        See function documentation for more details
#===============================================================================

def center_map(feature,**kwargs):
    """Center-mapping function used as an argument to :py:meth:`GenomeArray.add_from_bowtie`.
    A user-specified number of bases is optionally removed from each side of each
    read alignment, and the `N` remaining bases are each apportioned `1/N`
    of the read count, so that the entire read is counted once.
    
    Parameters
    ----------
    feature : |SegmentChain|
        Ungapped genomic alignment
        
    kwargs['value'] : float, optional
        Total value to divide over aligning positions (Default: `1.0`)
        
    kwargs['nibble'] : int, optional
        Positions to remove from each end before mapping (Default: `0`)
        
    kwargs['offset'] : int, optional
        Mapping offset, if any, from 5' end of read (Default: `0`)
    
    Returns
    -------
    list
        tuples of `(GenomicSegment,float value over segment)`
    """
    nibble = kwargs["nibble"]
    read_length = feature.length
    if read_length <= 2*nibble:
        warnings.warn("File contains read alignments shorter (%s nt) than `2*'nibble'` value of %s nt. Ignoring these." % (read_length,2*nibble),
                      DataWarning)
        return []

    span = feature.spanning_segment
    strand = span.strand
    chrom  = span.chrom
    offset = kwargs.get("offset",0)
    value  = float(kwargs.get("value",1.0))
    sign   = -1 if strand == "-" else 1
    start  = span.start + nibble + sign*offset
    end    = span.end   - nibble + sign*offset
    seg = GenomicSegment(chrom,
                         start,
                         end,
                         strand)
    frac = value/len(seg)
    return [(seg,frac)]

def five_prime_map(feature,**kwargs):
    """Fiveprime mapping function used as an argument to :py:meth:`GenomeArray.add_from_bowtie`.
    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
    
    Parameters
    ----------
    feature : |SegmentChain|
        Ungapped genomic alignment
        
    kwargs['value'] : float or int, optional
        Value to apportion (Default: `1`)
        
    kwargs['offset'] : int, optional
        Mapping offset, if any, from 5' toward 3' end of read
    
    Returns
    -------
    list
        tuples of `(GenomicSegment,float value over segment)`
    """
    offset  = kwargs.get("offset",0)
    read_length = feature.length
    if offset > read_length:
        warnings.warn("File contains read alignments shorter (%s nt) than offset (%s nt). Ignoring." % (read_length,offset),
                     DataWarning)
        return []

    value   = kwargs.get("value",1.0)
    span = feature.spanning_segment
    strand = span.strand
    if strand in ("+","."):
        start = span.start + offset
    else:
        start = span.end - 1 - offset
    seg = GenomicSegment(span.chrom,
                         start,
                         start+1,
                         strand)
    return [(seg,value)]

def three_prime_map(feature,**kwargs):
    """Threeprime mapping function used as an argument to :py:meth:`GenomeArray.add_from_bowtie`.
    Reads are mapped at a user-specified offset from the threeprime end of the alignment.
    
    Parameters
    ----------
    feature : |SegmentChain|
        Ungapped genomic alignment
        
    kwargs['value'] : float or int, optional
        Value to apportion (Default: `1`)
        
    kwargs['offset'] : int, optional
        Mapping offset, if any, from 3' toward 5' end of read.
    
    Returns
    -------
    list
        tuples of `(GenomicSegment,float value over segment)`
    """
    offset  = kwargs.get("offset",0)
    read_length = feature.length
    if offset > read_length:
        warnings.warn("File contains read alignments shorter (%s nt) than offset (%s nt). Ignoring." % (read_length,offset),
                      DataWarning)
        return []

    value   = kwargs.get("value",1.0)
    span = feature.spanning_segment
    strand = span.strand
    if strand in ("+","."):
        start = span.end - 1 - offset
    else:
        start = span.start + offset
    seg = GenomicSegment(span.chrom,
                         start,
                         start+1,
                         strand)
    return [(seg,value)]

def variable_five_prime_map(feature,**kwargs):
    """Fiveprime variable mapping function used as an argument to :py:meth:`GenomeArray.add_from_bowtie`.
    Reads are mapped at a user-specified offset from the fiveprime end of the alignment.
    The offset for a read of a given length is supplied in `kwargs[offset][readlen]`
    
    Parameters
    ----------
    feature : |SegmentChain|
        Ungapped genomic alignment
        
    kwargs['offset'] : dict, required
        Dictionary mapping read lengths to offsets

    kwargs['value'] : float or int, optional
        Value to apportion (Default: `1`)
    
    Returns
    -------
    list
        tuples of `(GenomicSegment,float value over segment)`
    """
    span = feature.spanning_segment
    strand = span.strand
    value   = kwargs.get("value",1.0)
    offset  = kwargs["offset"].get(len(span),kwargs["offset"].get("default",None))
    feature_length = feature.length
    if offset is None:
        warnings.warn("No offset for reads of length %s. Ignoring." % feature_length,
                      DataWarning)
        return []
    if offset >= feature_length:
        warnings.warn("Offset (%s nt) longer than read length %s. Ignoring" % (offset,feature_length))

    if strand in ("+","."):
        start = span.start + offset
    else:
        start = span.end - 1 - offset
    seg = GenomicSegment(span.chrom,
                         start,
                         start+1,
                         strand)
    return [(seg,value)]    

    

#===============================================================================
# GenomeArray classes
#===============================================================================

class AbstractGenomeArray(object):
    """Abstract base class for all |GenomeArray|-like objects"""
    
    def __init__(self,chr_lengths=None,strands=("+","-")):
        """
        Parameters
        ----------
        chr_lengths : dict
            Dictionary mapping chromosome names to lengths. Suppyling this in
            advance yields considerable speed improvements, as memory can be
            pre-allocated for each chromosomal position. If not provided, a
            minimum chromosome size will be guessed, and chromosomes re-sized
            as needed.
        
        strands : sequence, optional
            Strands for the |AbstractGenomeArray|. (Default: `("+","-")`)
        """
        self._chroms       = {}
        self._strands      = strands
        self._sum          = None
        self._normalize    = False

    def __str__(self):
        return repr(self)
    
    def __repr__(self):
        stmp = "<%s len=%s sum=%s chroms=" % (self.__class__.__name__, len(self), self.sum())
        stmp += ",".join(self.chroms())
        stmp += " strands=%s" % ",".join(self.strands())
        stmp += ">"
        return stmp

    def __contains__(self,chrom):
        """Return `True` if `chrom` is defined in the array, otherwise False
        
        Returns
        -------
        bool
        """
        return chrom in self._chroms

    def __len__(self):
        """Return size of GenomeArray in nucleotide positions (not base pairs).
        To obtain size in base pairs, divide length by two.
        
        Returns
        -------
        int
        """
        return len(self._strands)*sum(self.lengths().values())
    
    @abstractmethod
    def __getitem__(self,seg):
        """Retrieve array of counts from a region of interest. The values in
        the returned array are in 5' to 3' with respect to `seg` rather than
        the genome (i.e. are reversed for reverse-strand features).
        
        Parameters
        ----------
        seg : |GenomicSegment| or |SegmentChain|
            Region of interest in genome
        
        Returns
        -------
        :class:`numpy.ndarray`
            vector of numbers, each position corresponding to a position
            in `seg`, from 5' to 3' relative to `seg`
        
        See also
        --------
        SegmentChain.get_counts
            Fetch a spliced vector of data covering a |SegmentChain|
        """
        pass
    
    @abstractmethod    
    def reset_sum(self):
        """Reset sum to total mapped reads in the GenomeArray"""

    def set_sum(self,val):
        """Set sum used for normalization to an arbitrary value (e.g. from another dataset)
        
        Parameters
        ----------
        val : int or float
            a number
        """
        self._sum = val
        
    def sum(self):
        """Return the total number of aligned reads or the sum of the 
        quantitative data across all positions in the GenomeArray
        
        Notes
        -----
        The true (i.e. unnormalized) sum is always reported, even if 
        :meth:`set_normalize` is set to True
        
        Returns
        -------
        int or float
        """
        if self._sum is None:
            self.reset_sum()
        
        return self._sum

    def set_normalize(self,value=True):
        """Toggle normalization of reported values to reads per million mapped in the dataset
        
        Parameters
        ----------
        value : bool
            If `True`, all values fetched will be normalized to reads
            per million. If `False`, all values will not be normalized.
        """
        assert value in (True,False)
        self._normalize = value

    def chroms(self):
        """Return a list of chromosomes in the GenomeArray
        
        Returns
        -------
        list
            Chromosome names as strings
        """
        return self._chroms.keys()
    
    def strands(self):
        """Return a tuple of strands in the GenomeArray
        
        Returns
        -------
        tuple
            Chromosome strands as strings
        """
        return self._strands

    def lengths(self):
        """Return a dictionary mapping chromosome names to lengths.
        When two strands report different lengths for a chromosome, the
        max length is taken.
        
        Returns
        -------
        dict
            Dictionary mapping chromosome names to chromosome lengths
        """
        d_out = {}.fromkeys(self.keys())
        for key in d_out:
            d_out[key] = max([len(self._chroms[key][X]) for X in self.strands()])
        return d_out


class MutableAbstractGenomeArray(AbstractGenomeArray):
    """Abstract base class for |GenomeArray|-like objects whose values can be
    changed
    """
    @abstractmethod
    def __setitem__(self,seg,val):
        """Set values in |MutableAbstractGenomeArray| over a region of interest.
        
        Parameters
        ----------
        seg : |GenomicSegment| or |SegmentChain|
            Region of interest

        val : int, float, or :class:`numpy.ndarray`
            Scalar or vector of values to set in array over `seg`.
            If a vector, values should be ordered 5'-to-3' relative
            to `seg` rather than (i.e. for a reverse-strand feature,
            position 0 in the vector would correspond to seg.end)
            in the genome
        """
        pass


class BAMGenomeArray(AbstractGenomeArray):
    """A GenomeArray for :term:`read alignments` in `BAM`_ files.
    
    When a user requests :term:`read alignments` or :term:`counts` at
    a set of genomic positions, |BAMGenomeArray| uses a :term:`mapping rule`
    to determine which reads from the `BAM`_ file should be returned and how
    these reads should be mapped to positions in a vector. Reads
    are not mapped until requested, and mapping rules are changeable at
    runtime via :meth:`BAMGenomeArray.set_mapping`
    
    
    Attributes
    ----------
    map_fn : function
        :term:`mapping function` used to convert :term:`read alignments` to :term:`counts`
    
    bamfiles : list
        `BAM`_ files or :term:`read alignments`, as a list of :py:class:`pysam.AlignmentFile`

    
    Notes
    -----
    The alignment data in |BAMGenomeArray| objects is immutable. If you need
    to change the data in-place, the |BAMGenomeArray| can be converted to a
    |GenomeArray| or |SparseGenomeArray| via :meth:`~BAMGenomeArray.to_genome_array`
    """
    
    def __init__(self,bamfiles,mapping=CenterMapFactory()):
        """Create a |BAMGenomeArray|
        
        Parameters
        ----------
        bamfile : list
            An list of open :py:class:`pysam.AlignmentFile` s. Note: the
            corresponding `BAM`_ files must be sorted and indexed by `samtools`_. 
        
        mapping : func
            :term:`mapping function` that determines how each read alignment is mapped to a 
            count at a genomic position. Examples include mapping reads to their
            fiveprime ends, mapping reads to their threeprime ends, or mapping
            somewhere in the middle. Factories to produce such functions are provided.
            See references below. (Default: :func:`CenterMapFactory`)

        See Also
        --------
        FivePrimeMapFactory
            map reads to 5' ends, with or without applying an offset
        
        VariableFivePrimeMapFactory
            map reads to 5' ends, choosing an offset determined by read length
        
        ThreePrimeMapFactory
            map reads to 3' ends, with or without applying an offset
        
        CenterMapFactory
            map each read fractionally to every position in the read, optionally trimming positions from the ends first
        """
        self.bamfiles     = bamfiles
        self.map_fn       = mapping
        self._strands     = ("+","-")
        self._normalize   = False
        
        self._chr_lengths = {}
        for bamfile in self.bamfiles:
            for k,v in zip(bamfile.references,bamfile.lengths):
                self._chr_lengths[k] = max(self._chr_lengths.get(k,0),v)
                
        self._filters     = OrderedDict()
        self._update()

    def __del__(self):
        for bamfile in self.bamfiles:
            bamfile.close()

    def reset_sum(self):
        """Reset the sum to the total number of mapped reads in the |BAMGenomeArray|
        
        Notes
        -----
        Filters are not applied in this summation. See :meth:`BAMGenomeArray.set_sum`
        to manually set a sum, which e.g. could correspond to a total number of filtered
        reads.
        """
        self._sum = sum([X.mapped for X in self.bamfiles])
        
    def _update(self):
        """Updates mapping function to suit mapping rules
        """
        self.reset_sum()

    def add_filter(self,name,func):
        """Apply a function to filter reads retrieved from regions before mapping and counting
        
        Parameters
        ----------
        name : str
            A name for the filter. If not unique, will overwrite
            previous filter
        
        func : func
            Filter function. Function must take a
            :class:`pysam.AlignedSegment` as a single parameter, and return
            `True` if that read should be included in output, or `False` otherwise
        
        Notes
        -----
        In Python, `lambda` functions do NOT have their own scope! We strongly
        recomend defining filter functions using the ``def`` syntax to avoid
        namespace collisions.
        
        See also
        --------
        SizeFilterFactory
            generate filter functions that gate read alignments on size
        """
        self._filters[name] = func
    
    def remove_filter(self,name):
        """Remove a generic filter
        
        Parameters
        ----------
        name : str
            A name for the filter
        
        Returns
        -------
        func
            the removed filter function
        """
        retval = self._filters.pop(name)
        return retval
    
    def chroms(self):
        """Returns a list of chromosomes
        
        Returns
        -------
        list
            chromosome names as strings
        """
        return self._chr_lengths.keys()

    def lengths(self):
        """Returns a dictionary mapping chromosome names to lengths. 
        
        Returns
        -------
        dict
            mapping chromosome names to lengths
        """
        return self._chr_lengths
        
    def get_reads_and_counts(self,roi,roi_order=True):
        """Return :term:`read alignments` covering a |GenomicSegment|, and a
        count vector mapping reads to each positions in the |GenomicSegment|,
        following the rule specified by :meth:`~BAMGenomeArray.set_mapping`.
        Reads are strand-matched to `roi` by default. To obtain unstranded reads,
        set the value of `roi.strand` to `'.'`
        
        Parameters
        ----------
        roi : |GenomicSegment|
            Region of interest

        
        Returns
        -------
        list
            List of reads (as :class:`pysam.AlignedSegment`)
            covering region of interest

        numpy.ndarray
            Counts at each position of `roi`, under whatever mapping rule
            was set by :py:meth:`~BAMGenomeArray.set_mapping`


        Raises
        ------
        ValueError
            if bamfiles not sorted or not indexed
        """
        # fetch reads
        chrom  = roi.chrom
        strand = roi.strand
        start  = roi.start
        end    = roi.end
        if chrom not in self.chroms():
            return [], numpy.zeros(len(roi))

        reads = itertools.chain.from_iterable((X.fetch(reference=chrom,
                                               start=start,
                                               end=end,
                                               # until_eof=True, # this could speed things up. need to test/investigate
                                               ) for X in self.bamfiles))
            
        # filter by strand
        if strand == "+":
            reads = ifilter(lambda x: x.is_reverse is False,reads)
        elif strand == "-":
            reads = ifilter(lambda x: x.is_reverse is True,reads)
        
        # Pass through additional filters (e.g. size filters, if they have
        # been added)
        for my_filter in self._filters.values():
            reads = ifilter(my_filter,reads)
        
        # retrieve selected parts of regions
        reads,count_array = self.map_fn(list(reads),roi)
        
        # normalize to reads per million if normalization flag is set
        if self._normalize is True:
            count_array = count_array / float(self.sum()) * 1e6
        
        if roi_order == True and strand == "-":
            count_array = count_array[::-1]

        return reads, count_array

    def get_reads(self,roi):
        """Returns reads covering a |GenomicSegment|. Reads are strand-matched
        to `roi` by default, and are included or excluded depending upon the
        rule specified by :meth:`~BAMGenomeArray.set_mapping`. To obtain unstranded
        reads, set the value of `roi.strand` to `'.'`
        
        Parameters
        ----------
        roi : |GenomicSegment|
            Region of interest

        
        Returns
        -------
        list
            List of reads (as :class:`pysam.AlignedSegment`) covering region of interest,
            according to the rule set in :meth:`~BAMGenomeArray.set_mapping`


        Raises
        ------
        ValueError
            if bamfiles are not sorted or not indexed
        """
        reads, _ = self.get_reads_and_counts(roi)
        return reads

    def __getitem__(self,roi,roi_order=True): 
        """Retrieve array of counts from a region of interest, following
        the mapping rule set by :meth:`~BAMGenomeArray.set_mapping`.
        Values in the vector are ordered 5' to 3' relative to `roi`
        rather than the genome (i.e. are reversed for reverse-strand
        features).
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Region of interest in genome
        
        roi_order : bool, optional
            If `True` (default) return vector of values 5' to 3' 
            relative to vector rather than genome.

        Raises
        ------
        ValueError
            if bamfiles not sorted or not indexed

        Returns
        -------
        numpy.ndarray
            vector of numbers, each position corresponding to a position
            in `roi`, from 5' to 3' relative to `roi`
        
        See also
        --------
        SegmentChain.get_counts
            Fetch a spliced vector of data covering a |SegmentChain|
        """
        if isinstance(roi,SegmentChain):
            return roi.get_counts(self)

        _, count_array = self.get_reads_and_counts(roi,roi_order=roi_order)
        return count_array

    def get_mapping(self):
        """Return the docstring of the current mapping function
        """
        return self.map_fn.__doc__
    
    def set_mapping(self,mapping_function):
        """Change the mapping rule
        
        Parameters
        ----------
        mapping : func
            Function that determines how each read alignment is mapped to a 
            count at a position in the |BAMGenomeArray|. Examples include
            mapping reads to their fiveprime or threeprime ends, with or
            without offsets. Factories to produce such functions are provided.
            See references below.

            
        See Also
        --------
        FivePrimeMapFactory
            map reads to 5' ends, with or without applying an offset
        
        VariableFivePrimeMapFactory
            map reads to 5' ends, choosing an offset determined by read length
        
        ThreePrimeMapFactory
            map reads to 3' ends, with or without applying an offset
        
        CenterMapFactory
            map each read fractionally to every position in the read, optionally trimming positions from the ends first
        """
        self.map_fn = mapping_function
        self._update()
    
    def to_genome_array(self,array_type=None):
        """Converts |BAMGenomeArray| to a |GenomeArray| or |SparseGenomeArray|
        under the mapping rule set by :meth:`~BAMGenomeArray.set_mapping`

        Parameters
        ----------
        array_type : class
            Type of GenomeArray to return. |GenomeArray| or |SparseGenomeArray|.
            (Default: |GenomeArray|)
        
        Returns
        -------
        |GenomeArray| or |SparseGenomeArray|
        """
        if array_type is None:
            array_type = GenomeArray
            
        ga = array_type(chr_lengths=self.lengths(),strands=self.strands())
        for chrom in self.chroms():
            for strand in self.strands():
                seg = GenomicSegment(chrom,0,self.lengths()[chrom]-1,strand)
                ga[seg] = self[seg]

        return ga

    def to_variable_step(self,fh,trackname,strand,window_size=100000,**kwargs):
        """Write the contents of the |BAMGenomeArray| to a variableStep `Wiggle`_ file
        under the mapping rule set by :meth:`~BAMGenomeArray.set_mapping`.
        
        See the `Wiggle spec <http://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ for format details.
        
        Parameters
        ----------
        fh : file-like
            Filehandle to write to
        
        trackname : str
            Name of browser track
          
        strand : str
            Strand of |BAMGenomeArray| to export. `'+'`, `'-'`, or `'.'`
        
        window_size : int
            Size of chromosome/contig to process at a time.
            Larger values are faster but less memory-efficient
        
        **kwargs
            Any other key-value pairs to include in track definition line
        """
        assert strand in self.strands()
        fh.write("track type=wiggle_0 name=%s" % trackname)
        if kwargs is not None:
            for k,v in sorted(kwargs.items(),key = lambda x: x[0]):
                fh.write(" %s=%s" % (k,v))
        fh.write("\n")

        for chrom in sorted(self.chroms()):
            fh.write("variableStep chrom=%s span=1\n" % chrom)
            window_starts = xrange(0,self._chr_lengths[chrom],window_size)
            for i in range(len(window_starts)):
                my_start = window_starts[i]
                my_end   = window_starts[i+1] if i + 1 < len(window_starts) else int(self._chr_lengths[chrom])
                my_counts = self.__getitem__(GenomicSegment(chrom,
                                                            my_start,
                                                            my_end,
                                                            strand),
                                             roi_order=False)
                if my_counts.sum() > 0:
                    for idx in my_counts.nonzero()[0]:
                        genomic_x = my_start + idx
                        val = my_counts[idx]
                        fh.write("%s\t%s\n" % (genomic_x + 1,val))

            
        
    def to_bedgraph(self,fh,trackname,strand,window_size=100000,**kwargs):
        """Write the contents of the |BAMGenomeArray| to a `bedGraph`_ file
        under the mapping rule set by :meth:`~BAMGenomeArray.set_mapping`.
        
        See the `bedGraph spec <https://cgwb.nci.nih.gov/goldenPath/help/bedgraph.html>`_
        for format details.
            
        Parameters
        ----------
        fh : file-like
            Filehandle to write to
        
        trackname : str
            Name of browser track
          
        strand : str
            Strand of |BAMGenomeArray| to export. `'+'`, `'-'`, or `'.'`
        
        window_size : int
            Size of chromosome/contig to process at a time.
            Larger values are faster but less memory-efficient
        
        **kwargs
            Any other key-value pairs to include in track definition line
        """
        assert strand in self.strands()
        assert window_size > 0
        # write header
        fh.write("track type=bedGraph name=%s" % trackname)
        if kwargs is not None:
            for k,v in sorted(kwargs.items(),key = lambda x: x[0]):
                fh.write(" %s=%s" % (k,v))
        fh.write("\n")
        
        for chrom in sorted(self.chroms()):
            window_starts = xrange(0,self._chr_lengths[chrom],window_size)
            for i in range(len(window_starts)):
                my_start = window_starts[i]
                my_end   = window_starts[i+1] if i + 1 < len(window_starts) else int(self._chr_lengths[chrom])
                my_reads, my_counts = self.get_reads_and_counts(GenomicSegment(chrom,my_start,my_end,strand),roi_order=False)
                
                if len(my_reads) > 0:
                    genomic_start_x = window_starts[i]
                    last_val        = my_counts[0]

                    for x, val in enumerate(my_counts[1:]):
                        if val != last_val:
                            genomic_end_x = 1 + x + my_start
                            #write line: chrom chromStart chromEnd dataValue. 0-based half-open
                            if last_val > 0:
                                fh.write("%s\t%s\t%s\t%s\n" % (chrom,genomic_start_x,genomic_end_x,last_val))
                            
                            #update variables
                            last_val = val
                            genomic_start_x = genomic_end_x
                        else:
                            continue
                    # write out last values for window
                    if last_val > 0:
                        fh.write("%s\t%s\t%s\t%s\n" % (chrom,genomic_start_x,
                                                       my_end,
                                                       last_val))
        return


class GenomeArray(MutableAbstractGenomeArray):
    """Array-like data structure that maps numerical values (e.g. read :term:`counts`
    or conservation data, et c) to nucleotide positions in a genome.

    Data in the array may be manipulated manually and/or imported from `Wiggle`_
    or `BedGraph`_ files, as well as `bowtie`_ alignments.
    
    |GenomeArray| supports basic mathematical operations elementwise, where elements
    are nucleotide positions.
    """
    
    def __init__(self,chr_lengths=None,strands=("+","-"),
                 min_chr_size=MIN_CHR_SIZE):
        """Create a |GenomeArray|
        
        Parameters
        ----------
        chr_lengths : dict or `None`, optional
            Dictionary mapping chromosome names to lengths.
            Suppyling this in advance yields considerable
            speed improvements, as memory can be pre-allocated
            for each chromosomal position. If not provided,
            a minimum chromosome size will be guessed, and
            chromosomes re-sized as needed.
        
        min_chr_size : int
            If chr_lengths is not supplied, `min_chr_size` is 
            the default first guess of a chromosome size. If
            the genome has large chromosomes, it is much
            much better, speed-wise, to provide `chr_lengths`
            than to provide a guess here that is too small.
        
        strands : sequence
            Sequence of strand names for the |GenomeArray|. (Default: `('+','-')`)
        """
        self._chroms       = {}
        self._strands      = strands
        self.min_chr_size  = min_chr_size
        self._sum          = None
        self._normalize    = False
        if chr_lengths is not None:
            for chrom in chr_lengths.keys():
                self._chroms[chrom] = {}
                for strand in self._strands:
                    l = chr_lengths[chrom]
                    self._chroms[chrom][strand] = numpy.zeros(l)

    def reset_sum(self):
        """Reset the sum of the |GenomeArray| to the total number of mapped reads
        """
        self._sum = sum([X.sum() for X in self.iterchroms()])
        
    def _has_same_dimensions(self,other):
        """Return `True` if `self` and `other` have the chromosomes, strands, and chromosome lengths
        
        Parameters
        ----------
        other : |GenomeArray|
        
        
        Returns
        -------
        bool
        """
        assert set(self.keys()) == set(other.keys())
        assert self.strands() == other.strands()
        assert self.lengths() == other.lengths()

    def __eq__(self,other,tol=1e-10):
        """Test equality between `self` and `other`.
        
        To be equal, both |GenomeArrays| must have identical values at all
        nonzero positions in either array. Chromosome lengths need not be equal. 
        
        Parameters
        ----------
        other : |GenomeArray| or |SparseGenomeArray|
        
        
        Returns
        -------
        bool
        """
        snz = self.nonzero()
        onz = other.nonzero()
        for chrom in set(self.chroms()) | set(other.chroms()):
            for strand in set(self.strands()) | set(other.strands()):
                self_nonzero_vec  = snz.get(chrom,{K : numpy.array([]) for K in self.strands()}).get(strand,numpy.array([]))
                other_nonzero_vec = onz.get(chrom,{K : numpy.array([]) for K in other.strands()}).get(strand,numpy.array([]))

                if len(self_nonzero_vec) != len(other_nonzero_vec):
                    return False
                
                # indices of nonzero positions must be same
                if (self_nonzero_vec != other_nonzero_vec).any():
                    return False
                
                # values of nonzero positions must be same
                if len(self_nonzero_vec) > 0:
                    test_seg = GenomicSegment(chrom,self_nonzero_vec.min(),self_nonzero_vec.max(),strand)
                    if not (self[test_seg] - other[test_seg] <= tol).all():
                        return False

        return True
        
    def __getitem__(self,roi,roi_order=True):
        """Retrieve array of counts from a region of interest (`roi`)
        with values in vector ordered 5' to 3' relative to `roi`
        rather than genome (i.e. are reversed for reverse-strand
        features).
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Region of interest in genome
        
        roi_order : bool, optional
            If `True` (default) return vector of values 5' to 3' 
            relative to vector rather than genome.

        Returns
        -------
        :class:`numpy.ndarray`
            vector of numbers, each position corresponding to a position
            in `roi`, from 5' to 3' relative to `roi`
        
        See also
        --------
        SegmentChain.get_counts
            Fetch a spliced vector of data covering a |SegmentChain|
        """
        if isinstance(roi,SegmentChain):
            return roi.get_counts(self)
        
        chrom  = roi.chrom
        strand = roi.strand
        start  = roi.start
        end    = roi.end
        try:
            assert roi.end < len(self._chroms[chrom][strand])
        except AssertionError:
            my_len = len(self._chroms[chrom][strand])
            new_size = max(my_len + 10000,end + 10000)
            for my_strand in self.strands():
                new_strand = copy.deepcopy(self._chroms[roi.chrom][my_strand])
                new_strand.resize(new_size)
                self._chroms[chrom][my_strand] = new_strand
        except KeyError:
            assert strand in self.strands()
            if chrom not in self.keys():
                self._chroms[chrom] = {}
                for my_strand in self.strands():
                    self._chroms[chrom][my_strand] = numpy.zeros(self.min_chr_size)
        
        vals = self._chroms[chrom][strand][start:end]
        if self._normalize is True:
            vals = 1e6 * vals / self.sum()
            
        if roi_order == True and strand == "-":
            vals = vals[::-1]

        return vals
    
    def __setitem__(self,seg,val,roi_order=True):
        """Set values in the |GenomeArray| over a region of interest.
        
        If the `seg` is outside the bounds of the current |GenomeArray|,
        the |GenomeArray| will be automatically expanded to accomodate
        the genomic coordinates of `seg`.

        Parameters
        ----------
        seg : |GenomicSegment| or |SegmentChain|
            Region of interest

        val : int, float, or :class:`numpy.ndarray`
            Scalar or vector of values to set in array over `seg`.
            If a vector, values should be ordered 5'-to-3' relative
            to `seg` rather than (i.e. for a reverse-strand feature,
            position 0 in the vector would correspond to seg.end)
            in the genome

        roi_order : bool, optional
            If `True` (default) and `val` is a vector, values in `val`
            are assorted to go 5' to 3' relative to `seg` rather than
            genome
        """
        self._sum = None

        if isinstance(seg,SegmentChain):
            if isinstance(val,numpy.ndarray):
                if seg.spanning_segment.strand == "-":
                    val = val[::-1]

            length = seg.length
            x = 0
            for subseg in seg:
                if isinstance(val,numpy.ndarray):
                    subval = val[x:x+len(subseg)]
                else:
                    subval = val
                x += len(subseg)
                self.__setitem__(subseg,subval,roi_order=False)

            return

        old_normalize = self._normalize
        if old_normalize == True:
            warnings.warn("Temporarily turning off normalization during value set. It will be re-enabled automatically when complete.",DataWarning)

        self.set_normalize(False)


        chrom  = seg.chrom
        strand = seg.strand
        start  = seg.start
        end    = seg.end
        if strand == "-" and isinstance(val,numpy.ndarray) and roi_order == True:
            val = val[::-1]

        try:
            assert end < len(self._chroms[seg.chrom][seg.strand])
        except AssertionError:
            my_len = len(self._chroms[chrom][strand])
            new_size = max(my_len + 10000,end + 10000)
            for my_strand in self.strands():
                # this looks silly; but resize() can't work in-place iwth refs to array
                new_strand = copy.deepcopy(self._chroms[chrom][my_strand])
                new_strand.resize(new_size)
                self._chroms[chrom][my_strand]  = new_strand
        except KeyError:
            assert strand in self.strands()
            if chrom not in self.keys():
                self._chroms[chrom] = {}
                for my_strand in self.strands():
                    self._chroms[chrom][my_strand] = numpy.zeros(self.min_chr_size) 
            
        self._chroms[chrom][strand][start:end] = val
        self.set_normalize(old_normalize)

    def keys(self):
        """Return names of chromosomes in the GenomeArray"""
        return self.chroms()
    
    def iterchroms(self):
        """Return an iterator of each chromosome strand array
        
        Yields
        ------
        numpy.ndarray
            Positionwise counts on a chromosome strand
        """
        for chrom in self.keys():
            for strand in self.strands():
                yield self._chroms[chrom][strand]
    
    # no unit test for this
    def plot(self,chroms=None,strands=None,**plot_kw):
        """Create a plot of coverage along tracks or chromosomes
        
        Parameters
        ----------
        chroms : list
            Chromosomes to plot (default: all)
            
        strands : list
            Strands to plot (default: all)
            
        plot_kw
            A dictionary of keywords to pass to matplotlib
        
        
        Returns
        -------
        :py:class:`matplotlib.figure.Figure`
        """
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig_kw = { 'figsize' : (8,10) }
        colors = { "+" : "blue", "-" : "orange", "." : "darkgreen" }        
        
        if chroms is None:
            chroms = self.keys()
        
        if strands is None:
            strands = self.strands()
            
        fig_kw.update(plot_kw)
        f, axes = plt.subplots(nrows=len(chroms),**fig_kw)
        plt.subplots_adjust(hspace=0.3)
        for n,chrom in enumerate(sorted(chroms)):
            x = numpy.arange(0,self.lengths()[chrom])
            for strand in strands:
                multiplier = -1 if strand == "-" else 1
                axes[n].plot(x,self._chroms[chrom][strand]*multiplier,
                             label=strand,
                             color=colors[strand])
                axes[n].set_title(chrom)
                axes[n].set_ylabel("Counts")
                axes[n].set_xlabel("Position (nt)")
                axes[n].ticklabel_format(useOffset=1,style='plain',axis='x')
        
        return f
    
    def nonzero(self):
        """Return the indices of chromosomal positions with non-zero values
        at each chromosome/strand pair. Results are returned as a hierarchical
        dictionary mapping chromosome names, to a dictionary of strands,
        which in turn map to  to arrays of non-zero indices on that chromosome-strand.
        
        Returns
        -------
        dict
            `dict[chrom][strand]` = numpy.ndarray of indices
        """
        d_out = {}
        for key in self.keys():
            d_out[key] = {}
            for strand in self.strands():
                d_out[key][strand] = self._chroms[key][strand].nonzero()[0]
        return d_out
    
    def apply_operation(self,other,func,mode="same"):
        """Applies a binary operator to a copy of `self` and to `other` elementwise.
        `other` may be a scalar quantity or another |GenomeArray|. In both cases,
        a new |GenomeArray| is returned, and `self` is left unmodified.
        If :meth:`set_normalize` is set to `True`, it is disabled during the
        operation.
        
        Parameters
        ----------
        other : float, int, or |MutableAbstractGenomeArray|
            Second argument to `func`
        
        func : func
            Function to perform. This must take two arguments.
            a :py:class:`numpy.ndarray` (chromosome-strand) from the
            |GenomeArray| will be supplied as the first argument, and the
            corresponding chromosome-strand from `other` as the second.
            This operation will be applied over all chromosomes and strands.
            
            If mode is set to `all`, and a chromosome or strand is not
            present in one of self or other, zero will be supplied
            as the missing argument to func. func must handle this
            gracefully.
        
        mode : str, choice of `'same'`, `'all'`, or `'truncate'`
            Only relevant if `other` is a |MutableAbstractGenomeArray|
        
            If '`same`' each set of corresponding chromosomes must
            have the same dimensions in both parent |GenomeArray| s.
            In addition, both |GenomeArray| s must have the same
            chromosomes, and the same strands.

            If '`all`', corresponding chromosomes in the parent
            |GenomeArray| s will be resized to the dimensions of
            the larger chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; all chromosomes and strands will be
            included in output. 
                    
            If `'truncate'`,corresponding chromosomes in the parent
            |GenomeArray| s will be resized to the dimensions of
            the smaller chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; chromosomes or strands not common
            to both parents will be ignored.    
        
        Returns
        -------
        |GenomeArray|
            new |GenomeArray| after the operation is applied
        """
        new_array = GenomeArray.like(self)
        old_normalize = self._normalize
        if old_normalize == True:
            warnings.warn("Temporarily turning off normalization during value set. It will be re-enabled automatically when complete.",DataWarning)

        if type(other) == self.__class__:
            if mode == "same":       
                self._has_same_dimensions(other)
                for chrom in self.keys():
                    for strand in self.strands():
                        seg = GenomicSegment(chrom,0,len(self._chroms[chrom][strand]),strand)
                        new_array[seg] = func(self[seg],other[seg])               
            elif mode == "all":
                chroms    = {}.fromkeys(set(self.keys()) | set(other.keys()))
                strands   = set(self.strands()) | set(other.strands())
                new_array = GenomeArray(chroms,strands=strands)
                for chrom in chroms:
                    if chrom in self.keys() and chrom in other.keys():
                        for strand in strands:
                            if strand in self.strands() and strand in other.strands():
                                sc = copy.deepcopy(self._chroms[chrom][strand])
                                oc = copy.deepcopy(other._chroms[chrom][strand])
                                if len(sc) > len(oc):
                                    oc.resize(len(sc),refcheck=False)
                                elif len(oc) > len(sc):
                                    sc.resize(len(oc),refcheck=False)
                                new_array._chroms[chrom][strand] = func(sc,oc)
                            elif strand in self.strands():
                                new_array._chroms[chrom][strand] = func(copy.deepcopy(self._chroms[chrom][strand]),0)
                            else:
                                new_array._chroms[chrom][strand] = func(copy.deepcopy(other._chroms[chrom][strand]),0)
                    elif chrom in self.keys():
                        for strand in strands:
                            new_array._chroms[chrom][strand] = func(copy.deepcopy(self[chrom][strand]),0)
                    else:
                        for strand in strands:
                            new_array._chroms[chrom][strand] = func(copy.deepcopy(other[chrom][strand]),0)
            elif mode == "truncate":
                my_strands = set(self.strands()) & set(other.strands())
                my_chroms  = set(self.chroms())  & set(other.chroms())
                for chrom in my_chroms:
                    for strand in my_strands:
                        sc = copy.deepcopy(self._chroms[chrom][strand])
                        oc = copy.deepcopy(other._chroms[chrom][strand])
                        if len(sc) > len(oc):
                            sc.resize(len(oc),refcheck=False)
                        elif len(oc) > len(sc):
                            oc.resize(len(sc),refcheck=False)
                        new_array._chroms[chrom][strand] = func(sc,oc)
            else:
                raise ValueError("Mode not understood. Must be 'same', 'all', or 'truncate.'")
        else:
            for chrom in self.keys():
                for strand in self.strands():
                    new_array._chroms[chrom][strand] = func(self._chroms[chrom][strand],other)

        self.set_normalize(old_normalize)
        return new_array
        
    def __add__(self,other,mode="all"):
        """Add `other` to `self`, returning a new |GenomeArray|
        and leaving `self` unchanged. `other` may be a scalar quantity or
        another |GenomeArray|. In the latter case, addition is elementwise.
        
        Parameters
        ----------
        other : float, int, or |MutableAbstractGenomeArray|
        
        mode : str, choice of `'same'`, `'all'`, or `'truncate'`
            Only relevant if `other` is a |MutableAbstractGenomeArray|
        
            If '`same`' each set of corresponding chromosomes must
            have the same dimensions in both parent |GenomeArrays|.
            In addition, both |GenomeArrays| must have the same
            chromosomes, and the same strands.

            If '`all`', corresponding chromosomes in the parent
            |GenomeArrays| will be resized to the dimensions of
            the larger chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; all chromosomes and strands will be
            included in output. 
                    
            If `'truncate'`,corresponding chromosomes in the parent
            |GenomeArrays| will be resized to the dimensions of
            the smaller chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; chromosomes or strands not common
            to both parents will be ignored.     
        
        Returns
        -------
        |GenomeArray|
        """
        return self.apply_operation(other, operator.add, mode=mode)
        
    def __mul__(self,other,mode="all"):
        """Multiply `other` by `self`, returning a new |GenomeArray|
        and leaving `self` unchanged. `other` may be a scalar quantity or
        another |GenomeArray|. In the latter case, multiplication is elementwise.
        
        Parameters
        ----------
        other : float, int, or |MutableAbstractGenomeArray|
        
        mode : str, choice of `'same'`, `'all'`, or `'truncate'`
            Only relevant if `other` is a |MutableAbstractGenomeArray|
        
            If '`same`' each set of corresponding chromosomes must
            have the same dimensions in both parent |GenomeArrays|.
            In addition, both |GenomeArrays| s must have the same
            chromosomes, and the same strands.

            If '`all`', corresponding chromosomes in the parent
            |GenomeArrays| will be resized to the dimensions of
            the larger chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; all chromosomes and strands will be
            included in output. 
                    
            If `'truncate'`,corresponding chromosomes in the parent
            |GenomeArrays| will be resized to the dimensions of
            the smaller chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; chromosomes or strands not common
            to both parents will be ignored.                 
        
        Returns
        -------
        |GenomeArray|
        """
        return self.apply_operation(other, operator.mul, mode=mode)
        
    def __sub__(self,other,mode="all"):
        """Subtract `other` from `self`, returning a new |GenomeArray|
        and leaving `self` unchanged. `other` may be a scalar quantity or
        another |GenomeArray|. In the latter case, subtraction is elementwise.
        
        Parameters
        ----------
        other : float, int, or |MutableAbstractGenomeArray|
        
        mode : str, choice of `'same'`, `'all'`, or `'truncate'`
            Only relevant if `other` is a |MutableAbstractGenomeArray|
        
            If '`same`' each set of corresponding chromosomes must
            have the same dimensions in both parent |GenomeArrays|.
            In addition, both |GenomeArrays| s must have the same
            chromosomes, and the same strands.

            If '`all`', corresponding chromosomes in the parent
            |GenomeArrays| will be resized to the dimensions of
            the larger chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; all chromosomes and strands will be
            included in output. 
                    
            If `'truncate'`,corresponding chromosomes in the parent
            |GenomeArrays| will be resized to the dimensions of
            the smaller chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; chromosomes or strands not common
            to both parents will be ignored.                 
        
        Returns
        -------
        |GenomeArray|
        """
        return self.__add__(other*-1)
    
    def add_from_bowtie(self,fh,mapfunc,min_length=25,
                          max_length=numpy.inf,**trans_args):
        """Import alignment data in native `bowtie`_ format to the current GenomeArray
        
        Parameters
        ----------
        fh : file-like
            filehandle pointing to `bowtie`_ file
        
        mapfunc : func
            a :term:`mapping function` that transforms alignment
            coordinates into a sequence of (|GenomicSegment|,value) pairs

        min_length : int, optional
            minimum length for read to be counted
            (Default: 25)
        
        max_length : int or numpy.inf, optional
            maximum length for read to be counted
            (Default: infinity)
        
        trans_args : dict, optional
            Keyword arguments to pass to transformation
            function

            
        See also
        --------
        five_prime_map
            map reads to 5' ends, with or without applying an offset
        
        variable_five_prime_map
            map reads to 5' ends, choosing an offset determined by read length
        
        three_prime_map
            map reads to 3' ends, with or without applying an offset
        
        center_map
            map each read fractionally to every position in the read, optionally
            trimming positions from the ends first            
        """
        for feature in BowtieReader(fh):
            span_len = len(feature.spanning_segment)
            if span_len >= min_length and span_len <= max_length:
                tuples = mapfunc(feature,**trans_args)
                for seg, val in tuples:
                    self[seg] += val
        
        self._sum = None

    def add_from_wiggle(self,fh,strand):
        """Import data from a `Wiggle`_ or `bedGraph`_ file to current GenomeArray
        
        Parameters
        ----------
        fh : file-like
            filehandle pointing to wiggle file
        
        strand : str
            Strand to which data should be added. `'+'`, `'-'`, or `'.'`
        """
        assert strand in self.strands()
        for chrom,start,stop,val in WiggleReader(fh):
            seg = GenomicSegment(chrom,start,stop,strand)
            self[seg] += val

        self._sum = None
        
    def to_variable_step(self,fh,trackname,strand,**kwargs):
        """Export the contents of the GenomeArray to a variable step
        `Wiggle`_ file. For sparse data, `bedGraph`_ can be more efficient format.
        
        See the `Wiggle spec <http://genome.ucsc.edu/goldenpath/help/wiggle.html>`_
        for format details
        
        Parameters
        ----------
        fh : file-like
            Open filehandle to which data will be written
        
        trackname : str
            Name of browser track
        
        strand : str
            Strand to export. `'+'`, `'-'`, or `'.'`
        
        **kwargs
            Any other key-value pairs to include in track definition line
        """
        assert strand in self.strands()
        fh.write("track type=wiggle_0 name=%s" % trackname)
        if kwargs is not None:
            for k,v in sorted(kwargs.items(),key = lambda x: x[0]):
                fh.write(" %s=%s" % (k,v))
        fh.write("\n")
        nonzero = self.nonzero()
        for chrom in sorted(nonzero):
            fh.write("variableStep chrom=%s span=1\n" % chrom)
            indices = nonzero[chrom][strand]
            for idx in indices:
                val = self._chroms[chrom][strand][self._slicewrap(idx)]
                fh.write("%s\t%s\n" % (idx+1, val))
                
    def to_bedgraph(self,fh,trackname,strand,**kwargs):
        """Write the contents of the GenomeArray to a `bedGraph`_ file
        
        See the `bedGraph spec <https://cgwb.nci.nih.gov/goldenPath/help/bedgraph.html>`_
        for format details
            
        Parameters
        ----------
        fh : file-like
            Open filehandle to which data will be written
        
        trackname : str
            Name of browser track
        
        strand : str
            Strand to export. `'+'`, `'-'`, or `'.'`
        
        **kwargs
            Any other key-value pairs to include in track definition line
        """
        assert strand in self.strands()
        fh.write("track type=bedGraph name=%s" % trackname)
        if kwargs is not None:
            for k,v in sorted(kwargs.items(),key = lambda x: x[0]):
                fh.write(" %s=%s" % (k,v))
        fh.write("\n")
        nonzero = self.nonzero()
        for chrom in sorted(nonzero):
            if self._chroms[chrom][strand].sum() > 0:
                last_val = 0
                last_x   = 0
                nz = nonzero[chrom][strand]
                for x in range(nz.min(),nz.max()+1):
                    val = self._chroms[chrom][strand][self._slicewrap(x)]
                    if val != last_val:
                        #write line: chrom chromStart chromEnd dataValue
                        fh.write("%s\t%s\t%s\t%s\n" % (chrom,last_x,x,last_val))
                        #update variables
                        last_val = val
                        last_x = x
                    else:
                        continue
                # write last line
                fh.write("%s\t%s\t%s\t%s\n" % (chrom,last_x,x+1,last_val))
        
    @staticmethod
    def like(other):
        """Return a |GenomeArray| of same dimension as the input array
        
        Parameters
        ----------
        other : |GenomeArray|
        
        Returns
        -------
        GenomeArray
            empty |GenomeArray| of same size as `other`
        """
        return GenomeArray(other.lengths(),strands=other.strands())

    def _slicewrap(self,x):
        """Helper function to wrap coordinates for VariableStep/`bedGraph`_ export"""
        return x


class SparseGenomeArray(GenomeArray):
    """A memory-efficient sublcass of |GenomeArray| using sparse internal representation.
    Note, savings in memory may come at a cost in performance when repeatedly
    getting/setting values, compared to a |GenomeArray|.
    """
    def __init__(self,chr_lengths=None,strands=("+","-"),min_chr_size=MIN_CHR_SIZE):
        """Create a |SparseGenomeArray|

        Parameters
        ----------
        chr_lengths : dict, optional
            Dictionary mapping chromosome names to lengths.
            Suppyling this parameter yields considerable
            speed improvements, as memory can be pre-allocated
            for each chromosomal position. If not provided,
            a minimum chromosome size will be guessed, and
            chromosomes re-sized as needed (Default: `{}` )
        
        min_chr_size : int,optional
            If `chr_lengths` is not supplied, `min_chr_size` is 
            the default first guess of a chromosome size. If
            your genome has large chromosomes, it is much
            much better, speed-wise, to provide a chr_lengths
            dict than to provide a guess here that is too small.
            (Default: %s)
        
        strands : sequence
            Sequence of strand names for the |GenomeArray|. (Default: `('+','-')`)
        """ % MIN_CHR_SIZE
        self._chroms       = {}
        self._strands      = strands
        self._sum          = None
        self._normalize    = False
        self.min_chr_size = min_chr_size
        if chr_lengths is not None:
            for chrom in chr_lengths.keys():
                self._chroms[chrom] = {}
                for strand in self._strands:
                    l = chr_lengths[chrom]
                    self._chroms[chrom][strand] = scipy.sparse.dok_matrix((1,l))

    def lengths(self):
        """Return a dictionary mapping chromosome names to lengths. In the
        case where two strands report different lengths for a chromosome, the
        max length is taken.
        
        Returns
        -------
        dict
            mapping chromosome names to lengths
        """
        d_out = {}.fromkeys(self.keys())
        for key in d_out:
            d_out[key] = max([self._chroms[key][X].shape[1] for X in self.strands()])
        
        return d_out

    def __getitem__(self,roi,roi_order=True):
        """Retrieve array of counts from a region of interest (`roi`)
        with values in vector ordered 5' to 3' relative to `roi`
        rather than genome (i.e. are reversed for reverse-strand
        features).
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Region of interest in genome
        
        roi_order : bool, optional
            If `True` (default) return vector of values 5' to 3' 
            relative to `roi` rather than genome.

        Returns
        -------
        :class:`numpy.ndarray`
            vector of numbers, each position corresponding to a position
            in `roi`, from 5' to 3' relative to `roi`
        
        See also
        --------
        SegmentChain.get_counts
            Fetch a spliced vector of data covering a |SegmentChain|
        """
        if isinstance(roi,SegmentChain):
            return roi.get_counts(self)

        if roi.chrom not in self:
            self._chroms[roi.chrom] = { K : copy.deepcopy(scipy.sparse.dok_matrix((1,self.min_chr_size)))
                                       for K in self.strands()
                                      }
        if roi.end > self._chroms[roi.chrom][roi.strand].shape[1]:
            for strand in self.strands():
                self._chroms[roi.chrom][strand].resize((1,roi.end+10000))

        vals = self._chroms[roi.chrom][roi.strand][0,roi.start:roi.end]
        if self._normalize is True:
            vals = 1e6 * vals / self.sum()
            
        vals = vals.toarray().reshape(vals.shape[1],)
        if roi.strand == "-" and roi_order == True:
            vals = vals[::-1]

        return vals

    def __setitem__(self,seg,val,roi_order=True):
        """Set values in the |SparseGenomeArray| over a region of interest.
        
        If the `seg` is outside the bounds of the current |GenomeArray|,
        the |GenomeArray| will be automatically expanded to accomodate
        the genomic coordinates of `seg`.

        Parameters
        ----------
        seg : |GenomicSegment| or |SegmentChain|
            Region of interest

        val : int, float, or :class:`numpy.ndarray`
            Scalar or vector of values to set in array over `seg`.
            If a vector, values should be ordered 5'-to-3' relative
            to `seg` rather than (i.e. for a reverse-strand feature,
            position 0 in the vector would correspond to seg.end)
            in the genome

        roi_order : bool, optional
            If `True` (default) and `val` is a vector, values in `val`
            are assorted to go 5' to 3' relative to `seg` rather than
            genome
       """
        self._sum = None

        if isinstance(seg,SegmentChain):
            if isinstance(val,numpy.ndarray):
                if seg.spanning_segment.strand == "-":
                    val = val[::-1]

            length = seg.length
            x = 0
            for subseg in seg:
                if isinstance(val,numpy.ndarray):
                    subval = val[x:x+len(subseg)]
                else:
                    subval = val
                x += len(subseg)
                self.__setitem__(subseg,subval,roi_order=False)

            return

        old_normalize = self._normalize

        if isinstance(val,numpy.ndarray) and seg.strand == "-" and roi_order == True:
            val = val[::-1]

        if old_normalize == True:
            warnings.warn("Temporarily turning off normalization during value set. It will be re-enabled automatically when complete.",DataWarning)
        self.set_normalize(False)

        if seg.chrom not in self:
            self._chroms[seg.chrom] = { K : copy.deepcopy(scipy.sparse.dok_matrix((1,self.min_chr_size)))
                                       for K in self.strands()
                                      }
        if seg.end > self._chroms[seg.chrom][seg.strand].shape[1]:
            for strand in self.strands():
                self._chroms[seg.chrom][strand].resize((1,seg.end+10000))
        
        self._chroms[seg.chrom][seg.strand][0,seg.start:seg.end] = val
        self.set_normalize(old_normalize)

    def __mul__(self,other,mode=None):
        """Multiply `other` by `self`, returning a new |SparseGenomeArray|
        and leaving `self` unchanged. `other` may be a scalar quantity or
        another |SparseGenomeArray|. In the latter case, multiplication is elementwise.
        
        Parameters
        ----------
        other : float, int, or |MutableAbstractGenomeArray|
        
        mode : str, choice of `'same'`, `'all'`, or `'truncate'`
            Only relevant if `other` is a |MutableAbstractGenomeArray|
        
            If '`same`' each set of corresponding chromosomes must
            have the same dimensions in both parent |SparseGenomeArrays|.
            In addition, both |SparseGenomeArrays| must have the same
            chromosomes, and the same strands.

            If '`all`', corresponding chromosomes in the parent
            |SparseGenomeArrays| will be resized to the dimensions of
            the larger chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; all chromosomes and strands will be
            included in output. 
                    
            If `'truncate'`,corresponding chromosomes in the parent
            |SparseGenomeArrays| will be resized to the dimensions of
            the smaller chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; chromosomes or strands not common
            to both parents will be ignored.                 
        
        Returns
        -------
        |SparseGenomeArray|
        """
        if isinstance(other,GenomeArray):
            new_array = SparseGenomeArray.like(self)
            chroms    = set(self.keys()) & set(other.keys())
            strands   = set(self.strands()) & set(other.strands())
            for chrom in chroms:
                for strand in strands:
                    new_array._chroms[chrom][strand] = self._chroms[chrom][strand].multiply(other._chroms[chrom][strand])
                    
            return new_array
        else:
            return self.apply_operation(other,operator.mul,mode=mode)
        
    def apply_operation(self,other,func,mode=None):
        """Apply a binary operator to a copy of `self` and to `other` elementwise.
        `other` may be a scalar quantity or another |GenomeArray|. In both cases,
        a new |SparseGenomeArray| is returned, and `self` is left unmodified.
        If :meth:`set_normalize` is set to `True`, it is disabled during the
        operation.
        
        Parameters
        ----------
        other : float, int, or |MutableAbstractGenomeArray|
            Second argument to `func`
        
        func : func
            Function to perform. This must take two arguments.
            a :py:class:`numpy.ndarray` (chromosome-strand) from the
            |SparseGenomeArray| will be supplied as the first argument, and the
            corresponding chromosome-strand from `other` as the second.
            This operation will be applied over all chromosomes and strands.
            
            If mode is set to `all`, and a chromosome or strand is not
            present in one of self or other, zero will be supplied
            as the missing argument to func. func must handle this
            gracefully.
        
        mode : str, choice of `'same'`, `'all'`, or `'truncate'`
            Only relevant if `other` is a |MutableAbstractGenomeArray|
        
            If '`same`' each set of corresponding chromosomes must
            have the same dimensions in both parent |SparseGenomeArrays|.
            In addition, both |SparseGenomeArrays| must have the same
            chromosomes, and the same strands.

            If '`all`', corresponding chromosomes in the parent
            |SparseGenomeArrays| will be resized to the dimensions of
            the larger chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; all chromosomes and strands will be
            included in output. 
                    
            If `'truncate'`,corresponding chromosomes in the parent
            |SparseGenomeArrays| will be resized to the dimensions of
            the smaller chromosome before the operation is applied.
            Parents are not required to have the same chromosomes
            or strands; chromosomes or strands not common
            to both parents will be ignored.   
        
        Returns
        -------
        |SparseGenomeArray|
            new |SparseGenomeArray| after the operation is applied
        """
        new_array = SparseGenomeArray.like(self)

        old_normalize = self._normalize
        if old_normalize == True:
            warnings.warn("Temporarily turning off normalization during value set. It will be re-enabled automatically when complete.",DataWarning)
        self.set_normalize(False)

        if isinstance(other,GenomeArray):
            chroms    = {}.fromkeys(set(self.keys()) | set(other.keys()),10)
            strands   = set(self.strands()) | set(other.strands())
            new_array = SparseGenomeArray(chroms,strands=strands)
            for chrom in chroms:
                if chrom in self.keys() and chrom in other.keys():
                    for strand in strands:
                        if strand in self.strands() and strand in other.strands():
                            sc = self._chroms[chrom][strand]
                            oc = other._chroms[chrom][strand]
                            new_array._chroms[chrom][strand] = func(sc,oc)
                        elif strand in self.strands():
                            new_array._chroms[chrom][strand] = func(copy.deepcopy(self._chroms[chrom][strand]),0)
                        else:
                            new_array._chroms[chrom][strand] = func(copy.deepcopy(other._chroms[chrom][strand]),0)
                elif chrom in self.keys():
                    for strand in strands:
                        new_array._chroms[chrom][strand] = func(copy.deepcopy(self[chrom][strand]),0)
                else:
                    for strand in strands:
                        new_array._chroms[chrom][strand] = func(copy.deepcopy(other[chrom][strand]),0)
        else:
            for chrom in self.keys():
                for strand in self.strands():
                    new_array._chroms[chrom][strand] = func(self._chroms[chrom][strand],other)

        self.set_normalize(old_normalize)
        return new_array    
    
    def nonzero(self):
        """Return the indices of chromosomal positions with non-zero values
        at each chromosome/strand pair. Results are returned as a hierarchical
        dictionary mapping chromosome names, to a dictionary of strands,
        which in turn map to  to arrays of non-zero indices on that chromosome-strand.
        
        Returns
        -------
        dict
            `dict[chrom][strand]` = numpy.ndarray of indices
        """
        d_out = {}
        for key in self.keys():
            d_out[key] = {}
            for strand in self.strands():
                # need to sort because dok doens't guarantee sorted indices in nonzero()                
                d_out[key][strand] = numpy.array(sorted(self._chroms[key][strand].nonzero()[1]))
                
        return d_out

    def _slicewrap(self,x):
        """Helper function to wrap coordinates for VariableStep/`bedGraph`_ export"""
        return (0,x)

    @staticmethod
    def like(other):
        """Return a |SparseGenomeArray| of same dimension as the input array
        
        Parameters
        ----------
        other : |GenomeArray| or |SparseGenomeArray|
        
        Returns
        -------
        |SparseGenomeArray|
            of same size as `other`
        """
        return SparseGenomeArray(other.lengths(),strands=other.strands())





#===============================================================================
# Genome sizes for various assemblies
#===============================================================================

human = hg19 = {'chr1': 249250621,
 'chr10': 135534747,
 'chr11': 135006516,
 'chr11_gl000202_random': 40103,
 'chr12': 133851895,
 'chr13': 115169878,
 'chr14': 107349540,
 'chr15': 102531392,
 'chr16': 90354753,
 'chr17': 81195210,
 'chr17_ctg5_hap1': 1680828,
 'chr17_gl000203_random': 37498,
 'chr17_gl000204_random': 81310,
 'chr17_gl000205_random': 174588,
 'chr17_gl000206_random': 41001,
 'chr18': 78077248,
 'chr18_gl000207_random': 4262,
 'chr19': 59128983,
 'chr19_gl000208_random': 92689,
 'chr19_gl000209_random': 159169,
 'chr1_gl000191_random': 106433,
 'chr1_gl000192_random': 547496,
 'chr2': 243199373,
 'chr20': 63025520,
 'chr21': 48129895,
 'chr21_gl000210_random': 27682,
 'chr22': 51304566,
 'chr3': 198022430,
 'chr4': 191154276,
 'chr4_ctg9_hap1': 590426,
 'chr4_gl000193_random': 189789,
 'chr4_gl000194_random': 191469,
 'chr5': 180915260,
 'chr6': 171115067,
 'chr6_apd_hap1': 4622290,
 'chr6_cox_hap2': 4795371,
 'chr6_dbb_hap3': 4610396,
 'chr6_mann_hap4': 4683263,
 'chr6_mcf_hap5': 4833398,
 'chr6_qbl_hap6': 4611984,
 'chr6_ssto_hap7': 4928567,
 'chr7': 159138663,
 'chr7_gl000195_random': 182896,
 'chr8': 146364022,
 'chr8_gl000196_random': 38914,
 'chr8_gl000197_random': 37175,
 'chr9': 141213431,
 'chr9_gl000198_random': 90085,
 'chr9_gl000199_random': 169874,
 'chr9_gl000200_random': 187035,
 'chr9_gl000201_random': 36148,
 'chrM': 16571,
 'chrUn_gl000211': 166566,
 'chrUn_gl000212': 186858,
 'chrUn_gl000213': 164239,
 'chrUn_gl000214': 137718,
 'chrUn_gl000215': 172545,
 'chrUn_gl000216': 172294,
 'chrUn_gl000217': 172149,
 'chrUn_gl000218': 161147,
 'chrUn_gl000219': 179198,
 'chrUn_gl000220': 161802,
 'chrUn_gl000221': 155397,
 'chrUn_gl000222': 186861,
 'chrUn_gl000223': 180455,
 'chrUn_gl000224': 179693,
 'chrUn_gl000225': 211173,
 'chrUn_gl000226': 15008,
 'chrUn_gl000227': 128374,
 'chrUn_gl000228': 129120,
 'chrUn_gl000229': 19913,
 'chrUn_gl000230': 43691,
 'chrUn_gl000231': 27386,
 'chrUn_gl000232': 40652,
 'chrUn_gl000233': 45941,
 'chrUn_gl000234': 40531,
 'chrUn_gl000235': 34474,
 'chrUn_gl000236': 41934,
 'chrUn_gl000237': 45867,
 'chrUn_gl000238': 39939,
 'chrUn_gl000239': 33824,
 'chrUn_gl000240': 41933,
 'chrUn_gl000241': 42152,
 'chrUn_gl000242': 43523,
 'chrUn_gl000243': 43341,
 'chrUn_gl000244': 39929,
 'chrUn_gl000245': 36651,
 'chrUn_gl000246': 38154,
 'chrUn_gl000247': 36422,
 'chrUn_gl000248': 39786,
 'chrUn_gl000249': 38502,
 'chrX': 155270560,
 'chrY': 59373566}


drosophila = dm3 = {'2L': 23011544,
 '2LHet': 368872,
 '2R': 21146708,
 '2RHet': 3288761,
 '3L': 24543557,
 '3LHet': 2555491,
 '3R': 27905053,
 '3RHet': 2517507,
 '4': 1351857,
 'U': 10049037,
 'Uextra': 29004656,
 'X': 22422827,
 'XHet': 204112,
 'YHet': 347038,
 'dmel_mitochondrion_genome': 19517}

yeast = sgd2013 = saccer10 = sgdR64 = {'chrI': 230218,
 'chrII': 813184,
 'chrIII': 316620,
 'chrIV': 1531933,
 'chrIX': 439888,
 'chrV': 576874,
 'chrVI': 270161,
 'chrVII': 1090940,
 'chrVIII': 562643,
 'chrX': 745751,
 'chrXI': 666816,
 'chrXII': 1078177,
 'chrXIII': 924431,
 'chrXIV': 784333,
 'chrXV': 1091291,
 'chrXVI': 948066,
 'chrmt': 85779}



