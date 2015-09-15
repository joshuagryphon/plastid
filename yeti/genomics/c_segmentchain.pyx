#!/usr/bin/env python
"""This module defines object types that describe features in a genome or contig.


Important classes
-----------------
|GenomicSegment|
    A fundamental unit of a feature, similar to :py:class:`HTSeq.GenomicInterval`.
    |GenomicSegment| describes a single region of a genome, and is fully specified
    by a chromosome name, a start coordinate, an end coordinate, and a strand.
    
    |GenomicSegments| provide no feature annotation data, and are used
    primarily to construct |SegmentChains| or |Transcripts|, which do
    provide feature annotation data. |GenomicSegment| implements various methods
    to test equality to, overlap with, and containment of other |GenomicSegment|
    objects.

|SegmentChain|
    Base class for genomic features with rich annotation data. |SegmentChains|
    can contain zero or more |GenomicSegments|, and can therefore model
    discontinuous features -- such as multi-exon transcripts or gapped alignments --
    in addition to continuous features.  
    
    |SegmentChain| implements numerous convenience methods, e.g. for:

        - Converting coordinates between the genome and the spliced space of the
          |SegmentChain|
        
        - Fetching genomic sequence, read alignments, or count data over
          the |SegmentChain|, accounting for splicing of the segments and, for
          reverse-strand features, reverse-complementing of sequence

        - Slicing or fetching sub-regions of a |SegmentChain|
          
        - Testing for equality, inequality, overlap, containment, or coverage
          of other |SegmentChain| or |GenomicSegment| objects, in stranded 
          or unstranded manners
          
        - Exporting to `BED`_, `GTF2`_, or `GFF3`_ formats, for use with other
          software packages or within a genome browser.

    In addition, |SegmentChain| objects have attribute dictionaries that allow
    sotrage of arbitrary annotation information (e.g. gene IDs, GO terms, database
    cross-references, or miscellaneous notes).

|Transcript|
     Subclass of |SegmentChain| that adds convenience methods for fetching CDS,
     5' UTRs, and 3' UTRs, if the transcript is coding.


Examples
--------
Construct a |SegmentChain| from |GenomicSegments|::

    >>>
    >>>
    

Fetch a vector of spliced counts covering a |SegmentChain| from a |GenomeHash|::

    >>>
    >>>


Fetch a sub-region of a chain::

    >>>
    >>>


Find genomic coordinate of position 53 in a chain::

    >>>
    >>>


"""
__date__ = "2011-09-01"
__author__ = "joshua"
import re
import copy
import warnings
import array

import numpy
from numpy.ma import MaskedArray as MaskedArray
cimport numpy

from cpython cimport array

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from yeti.util.services.exceptions import DataWarning
from yeti.util.services.colors import get_str_from_rgb255, get_rgb255_from_str
from yeti.readers.gff_tokens import make_GFF3_tokens, \
                                    make_GTF2_tokens


from yeti.genomics.c_roitools cimport GenomicSegment, \
                                      positionlist_to_segments, \
                                      positions_to_segments, \
                                      strand_to_str, str_to_strand
from yeti.genomics.c_roitools import NullSegment
from yeti.genomics.c_common cimport ExBool, true, false, bool_exception, \
                                    Strand, forward_strand, reverse_strand,\
                                    unstranded, undef_strand


#===============================================================================
# constants
#===============================================================================

cdef hash_template = array.array('l',[]) 
cdef mask_template = array.array('i',[])

igvpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)")
segpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)\(([+-.])\)")
ivcpat = re.compile(r"([^:]*):([^(]+)\(([+-.])\)")

# compile time constants - __richcmp__ test
DEF LT  = 0
DEF LEQ = 1
DEF EQ  = 2
DEF NEQ = 3
DEF GT  = 4
DEF GEQ = 5

INT    = numpy.int
FLOAT  = numpy.float
DOUBLE = numpy.double
LONG   = numpy.long

ctypedef numpy.int_t    INT_t
ctypedef numpy.float_t  FLOAT_t
ctypedef numpy.double_t DOUBLE_t
ctypedef numpy.long_t   LONG_t


#===============================================================================
# io functions
#===============================================================================

def _get_attr(dtmp):
    ltmp=[]
    for k,v in sorted(dtmp.items()):
        ltmp.append("%s='%s'" % (k,v))
    
    return ",".join(ltmp)

def _format_segmentchain(segchain):
    """Formats a |SegmentChain| as a string, which when ``eval`` ed, should reconstruct the |SegmentChain|.
    Used for creating test datasets.
    
    Parameters
    ----------
    segchain : |SegmentChain|
    
    Returns
    -------
    str
    """
    iv_ltmp = []
    for iv in segchain:
        iv_ltmp.append("GenomicSegment('%s',%s,%s,'%s')" % (iv.chrom,iv.start,iv.end,iv.strand))
    
    stmp = "SegmentChain(%s,%s)" % (",".join(iv_ltmp),
                                    _get_attr(segchain.attr)
                                                          )
    return stmp

def _format_transcript(tx):
    """Formats a |Transcript| as a string, which when ``eval`` ed, should reconstruct the |Transcript|.
    Used for creating test datasets.
    
    Parameters
    ----------
    tx : |Transcript|
    
    Returns
    -------
    str
    """
    iv_ltmp = []
    for iv in tx:
        iv_ltmp.append("GenomicSegment('%s',%s,%s,'%s')" % (iv.chrom,iv.start,iv.end,iv.strand))
    
    stmp = "Transcript(%s,ID='%s',cds_genome_start=%s,cds_genome_end=%s)" % (",".join(iv_ltmp),
                                                          tx.get_name(),
                                                          tx.attr["cds_genome_start"],
                                                          tx.attr["cds_genome_end"]
                                                          )
    return stmp


#===============================================================================
# sorting functions
#===============================================================================

def sort_segmentchains_lexically(segchain):
    """Key function to sort a list of |SegmentChains| lexically by genomic position,
    by (in order of precedence): chromosome, start, end, strand, length (in nucleotides)
    and name.

    Parameters
    ----------
    feature : |SegmentChain|


    Returns
    -------
    str
        Chromosome name
        
    int
        Leftmost coordinate of `segchain`
        
    int
        Rightmost coordinate of `segchain`
        
    str
        Chromosome strand (`'+'`,`'-'`, or `'.'`)
    
    int
        Length (nt) of segment chain
    
    str
        Name of `segchain`
    
    dict
        Attributes of `segchain`
    """
    if isinstance(segchain,SegmentChain):
        length = segchain.get_length()
        
    return (segchain.spanning_segment.chrom,
            segchain.spanning_segment.start,
            segchain.spanning_segment.end,
            segchain.spanning_segment.strand,
            length,
            segchain.get_name(),
            segchain.attr
           )

def sort_segmentchains_by_name(segchain):
    """Key function to sort |SegmentChains| by name"""
    return segchain.get_name()

#===============================================================================
# higher-order classes that handle multi-feature structures, like transcripts
# or alignments
#===============================================================================

cdef class SegmentChain(object):
    """Base class for genomic features. |SegmentChains| can contain zero or more
    |GenomicSegments|, and therefore can model discontinuous, features -- such
    as multi-exon transcripts or gapped alignments -- in addition,
    to continuous features.
    
    Numerous convenience functions are supplied for:
    
      - converting between coordinates relative to the genome and relative
        to the internal coordinates of a spliced |SegmentChain|
        
      - fetching genomic sequence, read alignments, or count data, accounting
        for splicing of the segments, and, in the case of reverse-strand features,
        reverse-complementing
      
      - slicing or fetching sub-regions of a |SegmentChain|
      
      - testing equality, inequality, overlap, containment, coverage of, or
        sharing of segments with other |SegmentChain| or |GenomicSegment| objects
    
      - import/export to `BED`_, `PSL`_, `GTF2`_, and `GFF3`_ formats,
        for use in other software packages or in a genome browser.
    
    Intervals are sorted from lowest to greatest starting coordinate on their
    reference sequence, regardless of strand. Iteration over the SegmentChain
    will yield intervals from left-to-right in the genome.
    

    Attributes
    ----------
    spanning_segment : |GenomicSegment|
        A |GenomicSegment| spanning the endpoints of the |SegmentChain|

    strand : str
        The chromosome strand (`'+'`, `'-'`, or `'.'`)

    chrom : str
        Name of the chromosome on which the |SegmentChain| resides

    attr : dict
        Any miscellaneous attributes or annotation data


    See Also
    --------
	Transcript
        Transcript subclass, additionally providing richer `GTF2`_, `GFF3`_,
        and `BED`_ export, as well as methods for fetching coding regions
        and UTRs as subsegments
    """

    def __cinit__(self,*segments,**attr):
        self._segments        = []
        self._mask_segments   = []
        self._position_hash   = array.clone(hash_template,0,False)
        self._position_mask   = array.clone(mask_template,0,False)
        self._inverse_hash    = {}
        self.length = 0
        self.masked_length = 0
        self.spanning_segment = NullSegment

    def __init__(self,*segments,**attr):
        """Create an |SegmentChain| from zero or more |GenomicSegment| objects
        
        Example::
        
            >>> seg1 = GenomicSegment("chrI",2000,2500,"+")
            >>> seg2 = GenomicSegment("chrI",10000,11000,"+")
            >>> chain = SegmentChain(seg1,seg2,ID="example_chain",type="mRNA")
            
        
        Parameters
        ----------
        *segments : |GenomicSegment|
            0 or more GenomicSegments on the same strand
        
        **attr : dict
            Arbitrary attributes, including, for example:
        
            ====================    ============================================
            **Attribute**           **Description**
            ====================    ============================================
            ``type``                A feature type used for `GTF2`_/`GFF3`_ export
                                    of each interval in the |SegmentChain|. (Default: `'exon'`)
            
            ``ID``                  A unique ID for the |SegmentChain|.

            ``transcript_id``       A transcript ID used for `GTF2`_ export

            ``gene_id``             A gene ID used for `GTF2`_ export
            ====================    ============================================
        """
        self.attr = attr
        if "type" not in attr:
            self.attr["type"] = "exon"
        
        self.c_add_segments(segments)

    property length:
        def __get__(self):
            return self.length

    property masked_length:
        def __get__(self):
            return self.masked_length

    property spanning_segment:
        def __get__(self):
            return self.spanning_segment

    property chrom:
        """Chromosome the SegmentChain resides on"""
        def __get__(self):
            return self.spanning_segment.chrom

    property strand:
        """Strand of the SegmentChain"""
        def __get__(self):
            return strand_to_str(self.spanning_segment.c_strand)

    property c_strand:
        def __get__(self):
            return self.spanning_segment.c_strand

    property _segments:
        """List of |GenomicSegments| that comprise the |SegmentChain|"""
        def __get__(self):
            return self._segments

    property _mask_segments:
        """Trimmed |GenomicSegments| representing masked regions of the |SegmentChain|"""
        def __get__(self):
            return self._mask_segments

    property attr:
        def __get__(self):
            return self.attr
        def __set__(self,attr):
            self.attr = attr

    cpdef void _update(self):
        """Rebuilt position hashes and spanning segment, after adding positions
        """
        cdef:
            int num_segs = len(self)
            list segs
            GenomicSegment seg0
        
        self.sort()
        self._get_position_hash()
        self.length = sum([len(X) for X in self._segments])

        self._position_mask = array.clone(mask_template,0,False)
        self.add_masks()
       
        if num_segs == 0:
            self.spanning_segment = NullSegment 
        elif num_segs == 1:
            self.spanning_segment = self._segments[0]
        elif num_segs >1:
            segs = self._segments
            seg0 = segs[0]
            self.spanning_segment = GenomicSegment(seg0.chrom,
                                                   seg0.start,
                                                   segs[-1].end,
                                                   seg0.strand)

        else:
            raise RuntimeError("Segmentchain '%s' has negative intervals (%s)?" % (self.get_name(),num_segs))

    cdef void _get_position_hash(self):
        """Update position hash, mapping |SegmentChain| coordinates to genomic coordinates
        
        Returns
        -------
        long [:]
            Memoryview mapping each position in the |SegmentChain| to its 
            chromosomal coordinate
        """
        cdef list segments = self._segments
        cdef long length = sum([len(X) for X in segments])
        cdef long [:] my_hash = array.clone(hash_template,length,False)  
        cdef long x, c
        cdef GenomicSegment segment
        cdef dict ihash = {}

        c = 0
        for segment in self._segments:
            x = segment.start
            while x < segment.end:
                my_hash[c] = x
                ihash[x] = c
                c += 1
                x += 1
        
        self._position_hash = my_hash
        self._inverse_hash  = ihash

    cpdef void sort(self):
        """Sort component segments by ascending 5' chromosomal coordinate"""
        self._segments.sort()
        self._mask_segments.sort()

    def __repr__(self):
        sout = "<%s segments=%s" % (self.__class__.__name__, len(self))
        if len(self) > 0:
            sout += " bounds=%s:%s-%s(%s)" % (self[0].chrom,
                                              self[0].start,
                                              self[-1].end,
                                              self[0].strand)
            sout += " name=%s" % self.get_name()
        sout += ">"
        return sout

    def __str__(self):
        """String representation of |SegmentChain|. Inverse of :py:meth:`.from_str`
        Chains are represented as:
        
            `'chrom_name:segment1_start-segment1_end^segment2_start-segment2_end^...^(strand)'`
        
        Where all coordinates are zero-indexed and half-open. Masked segments
        are not saved in this representation; nor are attributes in `self.attr`
        """ 
        if len(self) > 0:
            ltmp = ["%s-%s" % (segment.start, segment.end) for segment in self]
            stmp = "^".join(ltmp)
            sout = "%s:%s(%s)" % (self.chrom,stmp,self.strand)
        else:
            sout = "na"
        return sout

    def __getitem__(self,index):
        """Fetch a |GenomicSegment| from the |SegmentChain|
        
        Parameters
        ----------
        index : int
            Index of interval to select, from left-to-right in genome
        
        Returns
        -------
        |GenomicSegment|
        """
        return self._segments[index]
           
    def __iter__(self):
        """Interation over each |GenomicSegment| in the |SegmentChain|,
        from left to right on the chromsome"""
        return iter(self._segments)
    
    def __next__(self):
        """Return next |GenomicSegment| in the |SegmentChain|,
        from left to right on the chromsome"""
        return next(self._segments)
    
    def next(self):
        """Return next |GenomicSegment| in the |SegmentChain|,
        from left to right on the chromsome"""
        return self.__next__()
    
    def __len__(self):
        """Return the number of |GenomicSegments| in the |SegmentChain|"""
        return len(self._segments)
    
    def __contains__(self, object other):
        """Tests whether |SegmentChain| contains another |SegmentChain|
        or |GenomicSegment|. Containment is defined for each type as follows:
        
        =====================   =======================================================
        **Type of `other`**     **True if**
        ---------------------   -------------------------------------------------------
        |SegmentChain|          If `self` and `other` both contain more than one segment,
                                all segment-segment junctions in `other` must be       
                                represented in `self`, in identical order, and all
                                positions covered by `other` must be present in `self`.
                                If `other` contains one segment, it must be fully
                                contained by one segment in `self`.                                  

        |GenomicSegment|        `other` must be completely contained        
                                within one of the segments in `self`                       
        =====================   =======================================================
        
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        
        Returns
        -------
        bool
        """
        cdef ExBool retval
        if isinstance(other,GenomicSegment):
            return self.__contains__(SegmentChain(other))
        elif isinstance(other,SegmentChain):
            retval = self.c_contains(other)
            if retval == bool_exception:
                raise RuntimeError("SegmentChain.unstranded_overlaps() errored on chain '%s' and object '%s'" % (self,other))
            return retval == true
        else:
            raise TypeError("The 'in'/containment operator is only defined for GenomicSegments and SegmentChains")
     
    cdef ExBool c_contains(self, SegmentChain other) except bool_exception:
        cdef:
            set selfpos, opos
            list myjuncs, ojuncs
            int i, mystart
            bint found

        if len(self) == 0 or len(other) == 0:
            return false
        elif self.chrom != other.chrom:
            return false
        elif self.c_strand != other.c_strand:
            return false
        elif other.length > self.length:
            return false
        elif len(other) == 1:
            for segment in self:
                if segment.contains(other[0]):
                    return true
            return false 
        else:
            selfpos = self.get_position_set()
            opos    = other.get_position_set()
            
            if opos & selfpos == opos:
                myjuncs = self.get_junctions()
                ojuncs  = other.get_junctions()
    
                found = False
                for i, myjunc in enumerate(myjuncs):
                    if ojuncs[0] == myjunc:
                        mystart = i
                        found   = True
                        break
            else:
                return false
            
            if found == True:
                return true if ojuncs == myjuncs[mystart:mystart+len(ojuncs)] else false
            else:
                return false

    def __richcmp__(self, object other, int cmpval):
        """Test whether `self` and `other` are equal or unequal. Equality is defined as
        identity of positions, chromosomes, and strands. Two |SegmentChain| with
        zero intervals, by convention, are not equal.
           
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        Returns
        -------
        bool
        """
        if isinstance(other,GenomicSegment):
            return self.__richcmp__(SegmentChain(GenomicSegment),cmpval) == true
        elif isinstance(other,SegmentChain):
            retval = self.c_richcmp(other,cmpval)
            if retval == bool_exception:
                raise RuntimeError("SegmentChain eq/ineq/et c comparison failed with chain '%s' on object '%'""" % (self,other))
            return retval == true

    cpdef ExBool c_richcmp(self, SegmentChain other, int cmpval) except bool_exception:
        cdef ExBool opval
        if cmpval == EQ:
            if len(self) == 0 or len(other) == 0:
                return false
            else:
                if self.chrom == other.chrom and\
                   self.c_strand == other.c_strand and\
                   self.get_position_set() == other.get_position_set():
                       return true
                return false
        elif cmpval == NEQ:
            opval = self.c_richcmp(other,EQ)
            if opval == bool_exception:
                pass
            elif opval == true:
                return false
            elif opval == false:
                return true
        else:
            raise ValueError("This operation is not defined for SegmentChains.")

    def shares_segments_with(self, object other):
        """Returns a list of |GenomicSegment| that are shared between `self` and `other`
           
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        Returns
        -------
        list
            List of |GenomicSegments| common to `self` and `other`
        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        if isinstance(other,GenomicSegment):
            tmp = SegmentChain(other)
            return self.shares_segments_with(tmp)
        elif isinstance(other,SegmentChain):
            return self.c_shares_segments_with(other)
        else:
            raise TypeError("SegmentChain.shares_segments_with() is defined only for GenomicSegments and SegmentChains. Found %s." % type(other))
 
    cdef list c_shares_segments_with(self, SegmentChain other):
        cdef:
            list shared
            GenomicSegment segment
        if self.chrom != other.chrom or self.c_strand != other.c_strand:
            return []
        else:
            shared = []
            for segment in other:
                if segment in self._segments:
                    shared.append(segment)
            return shared
   
    def  unstranded_overlaps(self, object other):
        """Return `True` if `self` and `other` share genomic positions
        on the same chromosome, regardless of their strands
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
         
        Returns
        -------
        bool
            `True` if `self` and `other` share genomic positions on the same
            chromosome, False otherwise. Strands of `self` and `other` need
            not match
            
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef ExBool retval
        if isinstance(other,GenomicSegment):
            return self.unstranded_overlaps(SegmentChain(other))
        elif isinstance(other,SegmentChain):
            retval = self.c_unstranded_overlaps(other)
            if retval == bool_exception:
                raise RuntimeError("SegmentChain.unstranded_overlaps() errored on chain '%s' and object '%s'" % (self,other))
            return retval == true
        else:
            raise TypeError("SegmentChain.unstranded_overlaps() is only defined for GenomicSegments and SegmentChains. Found %s." % type(other))

    cdef ExBool c_unstranded_overlaps(self, SegmentChain other) except bool_exception:
        cdef:
            set o_pos, my_pos
        if self.chrom != other.chrom:
            return false
        else:
            my_pos = self.get_position_set()
            o_pos  = other.get_position_set()
        if len(my_pos & o_pos) > 0:
            return true
        return false

    def overlaps(self, object other):
        """Return `True` if `self` and `other` share genomic positions on the same strand
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        Returns
        -------
        bool
            `True` if `self` and `other` share genomic positions on the same
            chromosome and strand; False otherwise.
        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef ExBool retval
        if isinstance(other,GenomicSegment):
            return self.overlaps(SegmentChain(other))
        elif isinstance(other,SegmentChain):
            retval = self.c_overlaps(other)
            if retval == bool_exception:
                raise RuntimeError("SegmentChain.overlaps() errored on chain '%s' and object '%s'" % (self,other))
            return retval == true
        else:
            raise TypeError("SegmentChain.overlaps() is defined only for GenomicSegments and SegmentChains. Found %s." % type(other))
 
    cdef ExBool c_overlaps(self, SegmentChain other) except bool_exception:
        cdef:
            Strand sstrand = self.c_strand
            Strand ostrand = other.c_strand
        if self.c_unstranded_overlaps(other) == true and self.c_strand & other.c_strand == other.c_strand:
            return true
        return false
    
    def antisense_overlaps(self, object other):
        """Returns `True` if `self` and `other` share genomic positions on opposite strands
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
         
        Returns
        -------
        bool
            `True` if `self` and `other` share genomic positions on the same
            chromosome but opposite strand; False otherwise.
                    
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef ExBool retval
        if isinstance(other,GenomicSegment):
            return self.antisense_overlaps(SegmentChain(other))
        elif isinstance(other,SegmentChain):
            retval = self.c_antisense_overlaps(other)
            if retval == bool_exception:
                raise RuntimeError("SegmentChain.antisense_overlaps() errored on chain '%s' and object '%s'" % (self,other))
            return retval == true
        else:
            raise TypeError("SegmentChain.antisense_overlaps() is only defined for GenomicSegments and SegmentChains. Found %s." % type(other))

    cdef ExBool c_antisense_overlaps(self, SegmentChain other) except bool_exception:
        cdef:
            Strand sstrand = self.c_strand
            Strand ostrand = other.c_strand
        if self.c_unstranded_overlaps(other) == true:
            if sstrand == unstranded or ostrand == unstranded or sstrand != ostrand:
                return true
        return false

    def covers(self, object other):
        """Return `True` if `self` and `other` share a chromosome and strand,
        and all genomic positions in `other` are present in `self`.
        By convention, zero-length |SegmentChains| are not covered by other
        chains.
        
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
         
        Returns
        -------
        bool
            `True` if `self` and `other` share a chromosome and strand, and all
            genomic positions in `other` are present in `self`. Otherwise `False`
        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef ExBool retval
        if isinstance(other,GenomicSegment):
            return self.covers(SegmentChain(other))
        elif isinstance(other,SegmentChain):
            retval = self.c_covers(other)
            if retval == bool_exception:
                raise RuntimeError("SegmentChain.covers() errored on chain '%s' and object '%s'" % (self,other))
            return retval == true
        else:
            raise TypeError("SegmentChain.covers() is only defined for GenomicSegments and SegmentChains. Found %s." % type(other))
 
    cdef ExBool c_covers(self, SegmentChain other) except bool_exception:
        """Return `True` if `self` and `other` share a chromosome and strand,
        and all genomic positions in `other` are present in `self`.
        By convention, zero-length |SegmentChains| are not covered by other
        chains.
        
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
         
        Returns
        -------
        bool
            `True` if `self` and `other` share a chromosome and strand, and all
            genomic positions in `other` are present in `self`. Otherwise `False`
        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        if len(self) == 0 or len(other) == 0:
            return false
        elif self.chrom  == other.chrom and\
             self.c_strand == other.c_strand and\
             other.get_position_set() & self.get_position_set() == other.get_position_set():
            return true
        return false
    
    cpdef SegmentChain get_antisense(self):
        """Returns an |SegmentChain| antisense to `self`, with empty `attr` dict.
        
        Returns
        -------
        SegmentChain
            |SegmentChain| antisense to `self`
        """
        cdef:
            Strand strand = self.spanning_segment.c_strand
            Strand new_strand
            list new_segments

        if strand == forward_strand:
            new_strand = reverse_strand
        elif strand == reverse_strand:
            new_strand = forward_strand
        elif strand == unstranded:
            new_strand = unstranded
        else:
            raise ValueError("Strand for SegmentChain '%s' is undefined and has no anti-sense." % self)

        new_segments = [GenomicSegment(X.chrom,X.start,X.end,strand_to_str(new_strand)) for X in self._segments]
        return SegmentChain(*new_segments)

    def get_position_list(self,copy = False):
        """Retrieve a sorted **end-inclusive** numpy array of genomic coordinates in this |SegmentChain|
         
        Parameters
        ----------
        copy : bool, optional
            If `False` (default), return a view of the |SegmentChain|'s
            internal position mapping. If `True`, return  a copy
 
        Returns
        -------
        :class:`numpy.ndarray`
            Ggenomic coordinates in `self`, as integers, in genomic order
        """
        return self.c_get_position_list(copy)

    cdef numpy.ndarray c_get_position_list(self, bint copy):
        """Retrieve a sorted end-inclusive list of genomic coordinates in this |SegmentChain|
        
        Parameters
        ----------
        copy : bool, optional
            If `False` (default), return a view of the |SegmentChain|'s
            internal position mapping. If `True`, return  a copy

        Returns
        -------
        :class:`numpy.ndarray`
            Ggenomic coordinates in `self`, as integers, in genomic order
        """
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._position_hash,dtype=long) 
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cpdef set get_position_set(self):
        """Retrieve an end-inclusive set of genomic coordinates included in this |SegmentChain|
        
        Returns
        -------
        set
            Set of genomic coordinates, as integers
        """
        return set(self.get_position_list())
   
    cpdef set get_masked_position_set(self):
        """Returns a set of genomic coordinates corresponding to positions in 
        `self` that **HAVE NOT** been masked using :meth:`SegmentChain.add_masks`

        Returns
        -------
        set
            Set of genomic coordinates, as integers
        """
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._position_hash,dtype=long)
        cdef numpy.ndarray[INT_t,ndim=1]  mask = numpy.asarray(self._position_mask,dtype=int)
        assert len(positions) == len(mask)
        return set(positions[mask == 0])
    
    def get_name(self):
        """Returns the name of this |SegmentChain|, first searching through
        `self.attr` for the keys `ID`, `Name`, and `name`. If no value is found
        for any of those keys, a name is generated using :meth:`SegmentChain.__str__`
        
        Returns
        -------
        str
            In order of preference, `ID` from `self.attr`, `Name` from
            `self.attr`, `name` from `self.attr` or ``str(self)`` 
        """
        name = self.attr.get("ID",
               self.attr.get("Name",
               self.attr.get("name",
                             str(self))))
        return name
    
    def get_gene(self):
        """Return name of gene associated with |SegmentChain|, if any, 
        by searching through `self.attr` for the keys `gene_id` and `Parent`.
        If one is not found, a generated gene name for the SegmentChain is 
        made from :py:meth:`get_name`.

        Returns
        -------
        str
            Returns in order of preference, `gene_id` from `self.attr`, 
            `Parent` from `self.attr` or ``'gene_%s' % self.get_name()``
        """
        gene = self.attr.get("gene_id",
               self.attr.get("Parent",
               "gene_%s" % self.get_name()))
        if isinstance(gene,list):
            gene = ",".join(sorted(gene))
            
        return gene
    
    def get_length(self):
        """Return total length, in nucleotides, of `self`
        
        Returns
        -------
        int
        """
        return self.length #_position_hash.shape[0] #sum([len(X) for X in self])

    def get_masked_length(self):
        """Return the total length, in nucleotides, of positions in `self`
        that have not been masked using :meth:`SegmentChain.add_masks`
        
        Returns
        -------
        int
        """
        return self.masked_length
        
    cdef tuple check_segments(self, tuple segments):
        """Check that all segments have same chromosome name and strand.
        If `self` already has segments in it, all |GenomicSegments|
        in `segments` must additionally have the same chromosme and
        strand as `self.`

        If so, return those as a tuple. If not, raise ValueError.

        Parameters
        ----------
        segments : tuple
            Tuple of |GenomicSegments|

        Returns
        -------
        str
            chromosome name

        str
            chromosome strand ('+', '-', or '.')

        Raises
        ------
        ValueError
            If segments in `segments` have differing chromosomes or strands
            from `self` or, if `self` is 0-length, from each other
        """
        cdef:
            GenomicSegment seg, seg0
            str my_chrom
            str my_strand

        if len(segments) > 0:
            seg0 = segments[0]
            if len(self) > 0:
                my_chrom  = self.chrom
                my_strand = self.strand
            else:
                my_chrom  = seg0.chrom
                my_strand = seg0.strand

            for seg in segments:
                if seg.chrom != my_chrom:
                    raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
                if seg.strand != my_strand:
                    raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))

            return (my_chrom,my_strand)
        else:
            return (self.chrom,self.strand)

    cdef c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            set positions
            GenomicSegment seg

        if len(segments) > 0:
            my_chrom, my_strand = self.check_segments(segments)
            # TODO: figure out a more efficient way?  Comparing segments to each other not ideal because quadratic
            # add new positions
            positions = self.get_position_set()
            for seg in segments:
                positions |= set(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positions_to_segments(my_chrom,my_strand,positions)
            self._update()

    def add_segments(self,*segments):
        """Add 1 or more |GenomicSegments| to the |SegmentChain|. If there are
        already segments in the chain, the incoming segments must be 
        on the same strand and chromosome as all others present.

        Parameters
        ----------
        segments : |GenomicSegment|
            One or more |GenomicSegment| to add to |SegmentChain|
        """
        cdef:
            str my_chrom, my_strand
            set positions
            GenomicSegment seg

        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
        
    def add_masks(self,*mask_segments):
        """Adds one or more |GenomicSegment| to the collection of masks.
        Masks will be trimmed to the positions of the |SegmentChain|
        during addition.

        Parameters
        ----------
        mask_segments : |GenomicSegment|
            One or more segments, in genomic coordinates, covering positions to
            exclude from return values of :meth:`get_masked_position_set`,
            :meth:`get_masked_counts`, or :meth:`get_masked_length`
		
        See also
        --------
        SegmentChain.get_masks
        SegmentChain.get_masks_as_segmentchain
        SegmentChain.reset_masks
        """
        cdef:
            str my_chrom, my_strand
            set positions = set()
            GenomicSegment seg
            int [:] pmask = array.clone(mask_template,self.length,True)
            long i, coord
            int tmpsum = 0
            dict ihash

        my_chrom, my_strand = self.check_segments(mask_segments)
           
        # add new positions to any existing masks
        for segment in list(mask_segments) + self._mask_segments:
            positions |= set(range(segment.start,segment.end))
         
        # trim away non-overlapping masks
        positions &= self.get_position_set()
        
        # regenerate list of ivs from positions, in case some were doubly-listed
        self._mask_segments = positions_to_segments(my_chrom,my_strand,positions)

        ihash = self._inverse_hash
        for i in positions:
            coord = ihash[i]
            pmask[coord] = 1
            tmpsum += 1

        self.masked_length = self.length - tmpsum
        self._position_mask = pmask
    
    def get_masks(self):
        """Return masked positions as a list of |GenomicSegments|
        
        Returns
        -------
        list
            list of |GenomicSegments| representing masked positions
        
        See also
        --------
        SegmentChain.get_masks_as_segmentchain
        
        SegmentChain.add_masks
        
        SegmentChain.reset_masks
        """
        return self._mask_segments
    
    def get_masks_as_segmentchain(self):
        """Return masked positions as a |SegmentChain|
        
        Returns
        -------
        |SegmentChain|
            Masked positions
        
        See also
        --------
        SegmentChain.get_masks
        
        SegmentChain.add_masks

        SegmentChain.reset_masks
        """        
        return SegmentChain(*self._mask_segments)
    
    def reset_masks(self):
        """Removes masks added by :py:meth:`add_masks`

        See also
        --------
        SegmentChain.add_masks
        """
        self._mask_segments = []
        self.add_masks()
            
    def get_junctions(self):
        """Returns a list of |GenomicSegments| representing spaces
        between the |GenomicSegments| in `self` In the case of a transcript,
        these would represent introns. In the case of an alignment, these
        would represent gaps in the query compared to the reference.
        
        Returns
        -------
        list
            List of |GenomicSegments| covering spaces between the intervals in `self`
            (e.g. introns in the case of a transcript, or gaps in the case of
            an alignment)
        """
        juncs = []
        for i in range(len(self)-1):
            seg1, seg2 = self[i], self[i+1]
            juncs.append(GenomicSegment(seg1.chrom,
                                        seg1.end,
                                        seg2.start,
                                        seg1.strand))
        return juncs
    
    def as_gff3(self, str feature_type=None, bint escape=True, list excludes=[]):
        """Format a length-1 |SegmentChain| as a line of `GFF3`_ output.
        
        Because `GFF3`_ files permit many schemas of parent-child hierarchy,
        and in order to reduce confusion and overhead, attempts to export
        a multi-interval |SegmentChain| will raise an :py:obj:`AttributeError`.
        
        Instead, users may export the individual features from which the
        multi-interval |SegmentChain| was constructed, or construct features
        for them, setting *ID*, *Parent*, and *type* attributes following
        their own conventions.

         
        Parameters
        ----------
        feature_type : str
            If not None, overrides the `type` attribute of `self.attr`
        
        escape : bool, optional
            Escape tokens in column 9 of `GFF3`_ output (Default: `True*)
        
        excludes : list, optional
            List of attribute key names to exclude from column 9
            (Default: `[]`)
        
        Returns
        -------
        str
            Line of `GFF3`_-formatted text
        
            
        Raises
        -----
        AttributeError
            if the |SegmentChain| has multiple intervals
            
        Notes
        -----
        Columns of `GFF3`_ are as follows
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if coding) 
                9     Attributes                       
            ======== =========

        For futher information, see
            - `GFF3 file format specification <http://www.sequenceontology.org/gff3.shtml>`_
            - `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`
            - `SO releases <http://sourceforge.net/projects/song/files/SO_Feature_Annotation/>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_            
        """
        cdef dict gff_attr
        cdef list always_excluded, ltmp

        if len(self) == 0: # empty SegmentChain
            return ""
        elif len(self) > 1:
            raise AttributeError("Attempted export of multi-interval %s" % self.__class__)
            
        gff_attr = copy.deepcopy(self.attr)
        feature_type = self.attr["type"] if feature_type is None else feature_type
        
        always_excluded = ["source",
                           "score",
                           "phase",
                           "cds_genome_start",
                           "cds_genome_end",
                           "thickstart",
                           "thickend",
                           "type",
                           "_bedx_column_order"]

        for segment in self:
            ltmp = self._get_8_gff_columns(segment,feature_type) +\
                   [make_GFF3_tokens(gff_attr,
                                     excludes=always_excluded+excludes,
                                     escape=True)]

        return "\t".join(ltmp) + "\n"
    
    def as_gtf(self, str feature_type=None, bint escape=True, list excludes=[]):
        """Format |SegmentChain| as a block of `GTF2`_ output.
        
        The `frame` or `phase` attribute (`GTF2`_ column 8) is valid only for `'CDS'`
        features, and, if not present in `self.attr`, is calculated assuming
        the |SegmentChain| contains the entire coding region. If the |SegmentChain|
        contains multiple intervals, the `frame` or `phase` attribute will
        *always* be recalculated.
        
        All attributes in `self.attr`, except those created upon import,
        will be propagated to all of the features that are generated.
        
        Parameters
        ----------
        feature_type : str
            If not None, overrides the "type" attribute of ``self.attr``
        
        escape : bool, optional
            Escape tokens in column 9 of GTF output (Default: True)
        
        excludes : list, optional
            List of attribute key names to exclude from column 8
            (Default: *[]*)
        
        Returns
        -------
        str
            Block of GTF2-formatted text
        
        
        Notes
        -----
        `gene_id` and `transcript_id` are required
            The `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_ requires
            that attributes `gene_id` and `transcript_id` be defined. If these
            are not present in `self.attr`, their values will be guessed 
            following the rules in :py:meth:`SegmentChain.get_gene` and 
            :py:meth:`SegmentChain.get_name`, respectively.
        
        Beware of attribute loss
            To save memory, only the attributes shared by all of the individual
            sub-features (e.g. exons) that were used to assemble this |Transcript|
            have been stored in `self.attr`. This means that upon re-export to `GTF2`_,
            these sub-features will be lacking any attributes that were specific
            to them individually. Formally, this is compliant with the 
            `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_, which states
            explicitly that only the attributes `gene_id` and `transcript_id`
            are supported.
            
        Columns of `GTF2`_ are as follows
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if coding) 
                9     Attributes. "gene_id" and "transcript_id" are required                        
            ======== =========
        
        For more info
            - `GTF2 file format specification <http://mblab.wustl.edu/GTF22.html>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_           
        """
        cdef:
            dict gtf_attr
            list ltmp1, ltmp2, always_excluded

        if len(self) == 0:
            return ""
        
        gtf_attr = copy.deepcopy(self.attr)
        gtf_attr["transcript_id"] = self.attr.get("transcript_id",self.get_name())
        gtf_attr["gene_id"]       = self.attr.get("gene_id",self.get_gene())
        feature_type = self.attr["type"] if feature_type is None else feature_type
        
        ltmp1 = []
        
        always_excluded = ["source",
                           "Parent",
                           "score",
                           "phase",
                           "cds_genome_start",
                           "cds_genome_end",
                           "thickstart",
                           "thickend",
                           "type",
                           "color",
                           "_bedx_column_order"]
        
        for segment in self:
            ltmp2 = self._get_8_gff_columns(segment,feature_type) +\
                   [make_GTF2_tokens(gtf_attr,
                                     excludes=always_excluded+excludes,
                                     escape=escape)]
            
            ltmp1.append("\t".join(ltmp2))            

        return "\n".join(ltmp1) + "\n"
    
    def _get_8_gff_columns(self, GenomicSegment segment, str feature_type):
        """Format columns 1-8 of GFF/GTF2/GFF3 files.
        
        Parameters
        ----------
        segment : |GenomicSegment|
            Segment to export
        
        feature_type : str
            Type of feature (for column 3 of output)
        
        
        Notes
        ------
        Columns of GFF files are as follows:
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if CDS) 
                9     Attributes. Formatting depends on flavor of GFF                      
            ======== =========        
        """
        cdef:
            str chrom  = segment.chrom
            str strand = segment.strand
            str phase  = "."
            long segment_start = segment.start
            long segment_end   = segment.end
            long new_segment_start
            list ltmp

        if feature_type == "CDS":
            # use phase/frame if known for length-1 features
            # called "phase" in GFF3 conventions; "frame" in GTF2
            if len(self) == 1 and ("phase" in self.attr or "frame" in self.attr):
                phase = str(self.attr.get("phase",self.attr.get("frame")))
            # otherwise calculate
            else:
                new_segment_start = self.get_segmentchain_coordinate(chrom,segment_start,strand,stranded=False)
                phase = str((3 - (new_segment_start % 3)) % 3)
        
        ltmp = [chrom,
                self.attr.get("source","."),
                feature_type,
                str(segment_start + 1),
                str(segment_end),
                str(self.attr.get("score",".")),
                strand,
                phase]
        
        return ltmp
    
    def as_bed(self, thickstart=None, thickend=None, as_int=True, color=None, extra_columns=None):
        """Format |SegmentChain| as a string of BED12[+X] output.
        
        If the |SegmentChain| was imported as a `BED`_ file with extra columns,
        these will be output in the same order, after the `BED`_ columns.

        Parameters
        ----------
        thickstart : int or `None`, optional
            If not `None`, overrides the genome coordinate that starts thick
            plotting in genome browser found in `self.attr['thickstart']`
    
        thickend : int or None, optional
            If not None, overrides the genome coordinate that stops
            thick plotting in genome browser found in `self.attr['thickend']`

        as_int : bool, optional
            Force `score` to integer (Default: `True`)
   
        color : str or None, optional
            Color represented as RGB hex string.
            If not none, overrides the color in `self.attr['color']`
    
        extra_columns : None or list, optional
            If `None`, and the |SegmentChain| was imported using the `extra_columns`
            keyword of :meth:`~yeti.genomics.roitools.SegmentChain.from_bed`,
            the |SegmentChain| will be exported in BED 12+X format, in which
            extra columns are in the same order as they were upon import. If no extra columns
            were present, the |SegmentChain| will be exported a aa BED12 line.

            If a list of attribute names, these attributes will be exported as
            extra columns in order, overriding whatever happened upon import. 
            If an attribute name is not in the `attr` dict of the |SegmentChain|,
            it will be exported with the default empty value "".

            If an empty list, no extra columns will be exported; the |SegmentChain|
            will be formatted as a BED12 line.


        Returns
        -------
        str 
            Line of BED12[+X]-formatted text


        Notes
        -----
        BED12 columns are as follows:
            ======== =========
            Column   Contains
            ======== =========
               1     Contig or chromosome
               2     Start of first block in feature (0-indexed)
               3     End of last block in feature (half-open)
               4     Feature name
               5     Feature score
               6     Strand
               7     thickstart (in chromosomal coordinates)
               8     thickend (in chromosomal coordinates)
               9     Feature color as RGB tuple
               10    Number of blocks in feature
               11    Block lengths
               12    Block starts, relative to start of first block
            ======== =========

        For more details
            See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        """
        cdef:
            list ltmp
            GenomicSegment span
            
        if len(self) > 0:
            score = self.attr.get("score",0)
            span = self.spanning_segment
            try:
                score = float(score)
                if as_int is True:
                    score = int(round(score))
            except ValueError:
                score = 0
            except TypeError:
                score = 0
            
            try:
                color = get_rgb255_from_str(self.attr.get("color","#000000")) if color is None else color
                color = str(color).strip("(").strip(")").replace(" ","")
            except ValueError:
                color = self.attr.get("color","0,0,0") if color is None else color
            
            thickstart = self.attr.get("thickstart",span.start) if thickstart is None else thickstart
            thickend   = self.attr.get("thickend",span.start)   if thickend   is None else thickend
            
            ltmp = [span.chrom,
                    span.start,
                    span.end,
                    self.get_name(),
                    score,
                    span.strand,
                    thickstart,
                    thickend,
                    color,
                    len(self),
                    ",".join([str(len(X)) for X in self]) + ",",
                    ",".join([str(X.start - self[0].start) for X in self]) + ","
                   ]            

            if extra_columns is None:
                extra_columns = self.attr.get("_bedx_column_order",[])

            if len(extra_columns) > 0:
                ltmp.extend([self.attr.get(X,"") for X in extra_columns])
            
            return "\t".join([str(X) for X in ltmp]) + "\n"
        else:
            # SegmentChain with no intervals
            return ""
    
    def as_psl(self):
        """Formats |SegmentChain| as `PSL`_ (blat) output.
        
        Notes
        -----
        This will raise an :py:class:`AttributeError` unless the following
        keys are present and defined in `self.attr`, corresponding to the
        columns of a `PSL`_ file:
        
            ======  ===================================
            Column  Key
            ======  ===================================
                1   ``match_length``
                2   ``mismatches``
                3   ``rep_matches``
                4   ``N``
                5   ``query_gap_count``
                6   ``query_gap_bases``
                7   ``target_gap_count``
                8   ``target_gap_bases``
                9   ``strand``
                10  ``query_name``
                11  ``query_length``
                12  ``query_start``
                13  ``query_end``
                14  ``target_name``
                15  ``target_length``
                16  ``target_start``
                17  ``target_end``
                19  ``q_starts`` : list of integers
                20  ``l_starts`` : list of integers
            ======  ===================================
        
        These keys are defined only if the |SegmentChain| was created by
        :py:meth:`SegmentChain.from_psl`, or if the user has defined them.
        
        See the `PSL spec <http://pombe.nci.nih.gov/genome/goldenPath/help/blatSpec.html>`_
        for more information.
        
        
        Returns
        -------
        str
            PSL-representation of BLAT alignment

        
        Raises
        ------
        AttributeError
            If not all of the attributes listed above are defined
        """
        cdef:
            dict attr = self.attr
            list ltmp
            str block_sizes, q_starts, t_starts
        try:
            ltmp = [
                attr["match_length"],
                attr["mismatches"],
                attr["rep_matches"],
                attr["N"],
                attr["query_gap_count"],
                attr["query_gap_bases"],
                attr["target_gap_count"],
                attr["target_gap_bases"],
                attr["strand"],
                attr["query_name"],
                attr["query_length"],
                attr["query_start"],
                attr["query_end"],
                attr["target_name"],
                attr["target_length"],
                attr["target_start"],
                attr["target_end"],
                len(self),
            ]
   
            block_sizes = ",".join([str(len(X)) for X in self]) + ","
            q_starts = ",".join([str(X) for X in self.attr["q_starts"]]) + ","
            t_starts = ",".join([str(X) for X in self.attr["t_starts"]]) + ","
            
            ltmp.append(block_sizes)
            ltmp.append(q_starts)
            ltmp.append(t_starts)
            return "\t".join(str(X) for X in ltmp) + "\n"
        except KeyError:
            raise AttributeError("SegmentChains only support PSL output if all PSL attributes are defined in self.attr: match_length, mismatches, rep_matches, N, query_gap_count, query_gap_bases, strand, query_length, query_start, query_end, target_name, target_length, target_start, target_end")
        
    def get_segmentchain_coordinate(self, str chrom, long genomic_x, str strand, bint stranded = True):
        """Finds the |SegmentChain| coordinate corresponding to a genomic position
        
        Parameters
        ----------
        chrom : str
            Chromosome name
            
        genomic_x : int
            coordinate, in genomic space
            
        strand : str
            Chromosome strand (`'+'`, `'-'`, or `'.'`)
            
        stranded : bool, optional
            If `True`, coordinates are given in stranded space
            (i.e. from 5' end of chain, as one might expect for a transcript).
            If `False`, coordinates are given from the left end of `self`,
            regardless of strand. (Default: `True`)
        
        
        Returns
        -------
        int
            Position in |SegmentChain|
            
        Raises
        ------
        KeyError
            if position outside bounds of |SegmentChain|
        """
        if chrom != self.chrom:
            raise ValueError("get_segmentchain_coordinate: query chromosome '%s' does not match chain '%s'" % chrom, self)
        if strand != self.strand:
            raise ValueError("get_segmentchain_coordinate: query strand '%s' does not match chain '%s'" % strand, self)

        return self.c_get_segmentchain_coordinate(genomic_x,stranded)

    cdef long c_get_segmentchain_coordinate(self, long genomic_x, bint stranded) except -1:
        """Finds the |SegmentChain| coordinate corresponding to a genomic position
        
        Parameters
        ----------
        genomic_x : int
            coordinate, in genomic space
            
        stranded : bool, optional
            If `True`, coordinates are given in stranded space
            (i.e. from 5' end of chain, as one might expect for a transcript).
            If `False`, coordinates are given from the left end of `self`,
            regardless of strand. (Default: `True`)
        
        
        Returns
        -------
        int
            Position in |SegmentChain|
            
        Raises
        ------
        KeyError
            if position outside bounds of |SegmentChain|
        """
        cdef dict ihash = self._inverse_hash

        if self.c_strand == reverse_strand and stranded == True:
            return self.length - ihash[genomic_x] - 1
        else:
            return ihash[genomic_x]

    def get_genomic_coordinate(self,x,stranded=True):
        """Finds genomic coordinate corresponding to position `x` in `self`
        
        Parameters
        ----------
        x : int
            position of interest, relative to |SegmentChain|
            
        stranded : bool, optional
            If `True`, `x` is assumed to be in stranded space (i.e. counted from
            5' end of chain, as one might expect for a transcript). If `False`,
            coordinates assumed to be counted the left end of the `self`,
            regardless of the strand of `self`. (Default: `True`)
        
                             
        Returns
        -------
        str 
            Chromosome name
        
        long
            Genomic cordinate corresponding to position `x`
        
        str
            Chromosome strand (`'+'`, `'-'`, or `'.'`)
        
        
        Raises
        ------
        IndexError
            if `x` is outside the bounds of the |SegmentChain|
        """
        cdef long new_x = self.c_get_genomic_coordinate(x, stranded)
        if new_x == -1:
            raise RuntimeError("Cannot fetch coordinate %s from SegmentChain %s (length %s)" % (x,self,self.length))

        return self.chrom, new_x, self.strand

    cdef long c_get_genomic_coordinate(self, long x, bint stranded) except -1:
        """Finds genomic coordinate corresponding to position `x` in `self`
        
        Parameters
        ----------
        x : int
            position of interest, relative to |SegmentChain|
            
        stranded : bool
            If `True`, `x` is assumed to be in stranded space (i.e. counted from
            5' end of chain, as one might expect for a transcript). If `False`,
            coordinates assumed to be counted the left end of the `self`,
            regardless of the strand of `self`.

        .. note::

           Unlike :meth:`SegmentChain.get_genomic_coordinate`, this only returns
           a `long`, as opposed to a tuple of (chrom_name, position, chrom_strand)
           

        Returns
        -------
        long
            Genomic coordinate corresponding to position `x`


        See also
        --------
        SegmentChain.get_genomic_coordinate
        """
        cdef:
            long [:] positions = self._position_hash
            long length = self.length
            long orig_index = x
            Strand strand = self.c_strand

        if strand == reverse_strand and stranded == True:
            x = length - x - 1

        if x < 0 or x >= length:
            raise IndexError("Position %s is outside bounds [0,%s) of SegmentChain '%s'" % (orig_index, length, self.get_name()))

        return positions[x]

    def get_subchain(self, long start, long end, bint stranded=True, **extra_attr):
        """Retrieves a sub-|SegmentChain| corresponding a range of positions
        specified in coordinates relative this |SegmentChain|. Attributes in
        `self.attr` are copied to the child SegmentChain, with the exception
        of `ID`, to which the suffix `'subchain'` is appended.
        
        Parameters
        ----------
        start : int
            position of interest in SegmentChain coordinates, 0-indexed
            
        end : int
            position of interest in SegmentChain coordinates, 0-indexed 
            and half-open
            
        stranded : bool, optional
            If `True`, `start` and `end` are assumed to be in stranded space (i.e. counted from
            5' end of chain, as one might expect for a transcript). If `False`,
            they assumed to be counted the left end of the `self`,
            regardless of the strand of `self`. (Default: `True`)

        extra_attr : keyword arguments
            Values that will be included in the subchain's `attr` dict.
            These can be used to overwrite values already present.
                          
        Returns
        -------
        |SegmentChain|
            covering parent chain positions `start` to `end` of `self`
        
        
        Raises
        ------
        IndexError
            if `start` or `end` is outside the bounds of the |SegmentChain|

        TypeError
            if `start` or `end` is None
        """
        cdef SegmentChain chain = self.c_get_subchain(start,end,stranded)
        cdef dict old_attr = copy.deepcopy(self.attr)
        old_attr.update(extra_attr)
        chain.attr = old_attr
        chain.attr["ID"] = "%s_subchain" % self.get_name()
        return chain

    cdef SegmentChain c_get_subchain(self, long start, long end, bint stranded):
        """Similar to :meth:`SegmentChain.get_subchain` but does not copy `attr` dict
        
        Parameters
        ----------
        start : long
            position of interest in SegmentChain coordinates, 0-indexed
            
        end : long
            position of interest in SegmentChain coordinates, 0-indexed 
            and half-open
            
        stranded : bool
            If `True`, `start` and `end` are assumed to be in stranded space (i.e. counted from
            5' end of chain, as one might expect for a transcript). If `False`,
            they assumed to be counted the left end of the `self`,
            regardless of the strand of `self`.

        Returns
        -------
        |SegmentChain|
            covering parent chain positions `start` to `end`

        See also
        --------
        SegmentChain.get_subchain
        """
        cdef:
             long [:] positions = self._position_hash
             long length = self.length
             long tmp
             list segs

        if start is None:
            raise TypeError('start coordinate supplied is None. Expected int')
        elif end is None:
            raise TypeError('end coordinate supplied is None. Expected int')

        if stranded == True and self.c_strand == reverse_strand:
            tmp = end
            end   = length - start 
            start = length - tmp 
            
        positions = positions[start:end]
        segs = positions_to_segments(self.chrom,self.strand,positions)

        return SegmentChain(*segs)

    def get_counts(self,ga,stranded=True):
        """Return list of counts or values at each position in `self`
        
        Parameters
        ----------
        ga : non-abstract subclass of |AbstractGenomeArray|
            GenomeArray from which to fetch counts
            
        stranded : bool, optional
            If `True` and the SegmentChain is on the minus strand,
            count order will be reversed relative to genome so that the
            array positions march from the 5' to 3' end of the chain.
            (Default: `True`)
            
            
        Returns
        -------
        numpy.ndarray
            Array of counts from `gnd` covering `self`
        """
        cdef:
            long c, cend
            GenomicSegment seg
            numpy.ndarray[DOUBLE_t,ndim=1] count_array = numpy.zeros(self.length,dtype=DOUBLE)

        if len(self) == 0:
            warnings.warn("%s is a zero-length SegmentChain. Returning 0-length count vector." % self.get_name(),DataWarning)

        c = 0
        for seg in self:
            cend = c + len(seg)
            count_array[c:cend] = ga.__getitem__(seg,roi_order=False)
            c = cend

        if self.c_strand == reverse_strand and stranded is True:
            count_array = count_array[::-1]
            
        return count_array
    
    def get_masked_counts(self,ga,stranded=True,copy=False):
        """Return counts covering `self` in dataset `gnd` as a masked array, in transcript 
        coordinates. Positions masked by :py:meth:`SegmentChain.add_mask` 
        will be masked in the array
        
        Parameters
        ----------
        gnd : non-abstract subclass of |AbstractGenomeArray|
            GenomeArray from which to fetch counts
            
        stranded : bool, optional
            If true and the |SegmentChain| is on the minus strand,
            count order will be reversed relative to genome so that the
            array positions march from the 5' to 3' end of the chain.
            (Default: `True`)

        copy : bool, optional
            If `False` (defualt) returns a view of the data; so changing
            values in the view changes the values in the |GenomeArray|
            if it is mutable. If `True`, a copy is returned instead.
            
            
        Returns
        -------
		:py:class:`numpy.ma.masked_array`
        """
        cdef:
            numpy.ndarray[DOUBLE_t,ndim=1] counts = self.get_counts(ga)
            numpy.ndarray[INT_t,ndim=1] mask = numpy.asarray(self._position_mask,dtype=int)

        return MaskedArray(counts,mask=mask.astype(bool),copy=copy)
        
#        valid_positions = [self.get_segmentchain_coordinate(self.spanning_segment.chrom,X,
#                                                            self.spanning_segment.strand) 
#                           for X in self.get_masked_position_set()]
#        for x in valid_positions:
#            atmp.mask[x] = False
        
#        return atmp
    
    def get_sequence(self,genome,stranded=True):
        """Return spliced genomic sequence of |SegmentChain| as a string
        
        Parameters
        ----------
        genome : dict or :class:`twobitreader.TwoBitFile`
            Dictionary mapping chromosome names to sequences.
            Sequences may be strings, string-like, or :py:class:`Bio.Seq.SeqRecord` objects
       
        stranded : bool
            If `True` and the |SegmentChain| is on the minus strand,
            sequence will be reverse-complemented (Default: True)
            
            
        Returns
        -------
        str
            Nucleotide sequence of the |SegmentChain| extracted from `genome`
        """
        if len(self) == 0:
            warnings.warn("%s is a zero-length SegmentChain. Returning empty sequence." % self.get_name(),DataWarning)
            return ""

        else:
            chromseq = genome[self.spanning_segment.chrom]
            ltmp = [chromseq[X.start:X.end] for X in self]
            stmp = "".join([str(X.seq) if isinstance(X,SeqRecord) else X for X in ltmp])

            if self.strand == "-"  and stranded is True:
                stmp = str(Seq(stmp,generic_dna).reverse_complement())
            
        return stmp
    
    def get_fasta(self,genome,stranded=True):
        """Formats sequence of SegmentChain as FASTA output
        
        Parameters
        ----------
        genome : dict or :class:`twobitreader.TwoBitFile`
            Dictionary mapping chromosome names to sequences.
            Sequences may be strings, string-like, or :py:class:`Bio.Seq.SeqRecord` objects
       
        stranded : bool
            If `True` and the |SegmentChain| is on the minus strand,
            sequence will be reverse-complemented (Default: True)

            
        Returns
        -------
        str
            FASTA-formatted seuqence of |SegmentChain| extracted from `genome`
        """
        return ">%s\n%s\n" % (self.get_name(),self.get_sequence(genome,stranded=stranded))

    @staticmethod
    def from_str(str inp):
        """Create a |SegmentChain| from a string formatted by :py:meth:`SegmentChain.str`:
           
            `chrom:start-end^start-end(strand)`
           
        where '^' indicates a splice junction between regions specified
        by `start` and `end` and `strand` is '+', '-', or '.'. Coordinates are
        0-indexed and half-open.


        Parameters
        ----------
        inp : str
			String formatted in manner of :py:meth:`SegmentChain.str`
          
          
        Returns
        -------
        |SegmentChain|
        """
        cdef:
            str chrom, middle, strand, sstart, send, piece
            list segs
            long start, end

        if inp in ("na","nan","None:(None)","None","none",None) or isinstance(inp,float) and numpy.isnan(inp):
            return SegmentChain()
        else:
            chrom,middle,strand = ivcpat.search(inp).groups()
            segs = []
            for piece in middle.split("^"):
                sstart,send = piece.split("-")
                start = long(sstart)
                end = long(send)
                segs.append(GenomicSegment(chrom,start,end,strand))
            return SegmentChain(*segs)
        
    @staticmethod
    def from_bed(str line, extra_columns=0):
        """Create a |SegmentChain| from a line from a `BED`_ file.
        The `BED`_ line may contain 4 to 12 columns, per the specification.
        These will be auto-detected and parsed appropriately.
        
        See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        for more details.

        Parameters
        ----------
        line
            Line from a `BED`_ file, containing 4 or more columns

        extra_columns: int or list optional
            Extra, non-BED columns in :term:`BED X+Y`_ format file corresponding to feature
            attributes. This is common in `ENCODE`_-specific `BED`_ variants.
            
            if `extra-columns` is:
            
              - an :class:`int`: it is taken to be the
                number of attribute columns. Attributes will be stored in
                the `attr` dictionary of the |SegmentChain|, under names like
                `custom0`, `custom1`, ... , `customN`.

              - a :class:`list` of :class:`str`, it is taken to be the names
                of the attribute columns, in order, from left to right in the file.
                In this case, attributes in extra columns will be stored under
                there respective names in the `attr` dict.

              - a :class:`list` of :class:`tuple`, each tuple is taken
                to be a pair of `(attribute_name, formatter_func)`. In this case,
                the value of `attribute_name` in the `attr` dict of the |SegmentChain|
                will be set to `formatter_func(column_value)`.
            
            (Default: 0)

        Returns
        -------
        |SegmentChain|
        """
        cdef:
            int num_bed_columns, num_extra_columns
            dict attr
            list frags = []
            list items = line.strip("\n").split("\t")

            #list column_formatters, frag_sizes, frag_offsets
            #set types
            #int num_extra_columns, i
            #long chrom_start, chrom_end
            #str chrom, default_id, k
            #dict base_columns, attr
        
        if isinstance(extra_columns,int):
            if extra_columns < 0:
                raise ValueError("Cannot make SegmentChain from BED input: if an integer, extra_columns must be non-negative.")
            num_extra_columns = extra_columns
            column_formatters = [("custom%s" % X,str) for X in range(extra_columns)]
        elif isinstance(extra_columns,list):
            num_extra_columns = len(extra_columns)
            types = set([type(X) for X in extra_columns])
            if len(types) > 1:
                raise ValueError("List of `extra_columns` contains mixed types. Cannot parse.")
            elif str in types:
                column_formatters = [(X,str) for X in extra_columns]
            elif tuple in types:
                if all([len(X) == 2 for X in extra_columns]) == False:
                    raise ValueError("Cannot make SegmentChain from BED input: if a list, extra_columns must be a list of tuples of (column_name,formatter_func)")
                column_formatters = extra_columns
        else:
            raise TypeError("Cannot make SegmentChain from BED input: extra_columns must be an int or list. Got a %s" % type(extra_columns))
            
        num_bed_columns = len(items) - num_extra_columns
        if num_bed_columns < 3:
            raise ValueError("BED format requires at least 3 columns. Found only %s." % num_bed_columns)
        
        chrom         = items[0]
        chrom_start   = long(items[1])
        chrom_end     = long(items[2])
        strand = "." if num_bed_columns < 6 else items[5]
    
        default_id  = "%s:%s-%s(%s)" % (chrom,chrom_start,chrom_end,strand)
    
        # dict mapping optional bed column to tuple of (Name,default value)
        # these values are used if any optional columns 4-12 are ommited
        bed_columns = { 3 :  ("ID",         default_id,    str),
                        4 :  ("score",      numpy.nan,     float),
                        #5 :  ("strand",    ".", strand),
                        6 :  ("thickstart", None,          int),
                        7 :  ("thickend",   None,          int),
                        8 :  ("color",      "0,0,0",       str),
                        9 :  ("blocks",     "1",             int),
                        10 : ("blocksizes", str(chrom_end - chrom_start),str),
                        11 : ("blockstarts","0",             str),
                      }
    
        # set attr defaults in case we're dealing with BED4-BED9 format
        attr = { KEY : DEFAULT for KEY,DEFAULT,_ in bed_columns.values() }
    
        # populate attr with real values from BED columns that are present
        for i, tup in sorted(bed_columns.items()):
            if num_bed_columns > i:
                key     = tup[0]
                default = tup[1]
                func    = tup[2]
                try:
                    attr[key] = func(items[i])
                except ValueError:
                    attr[key] = default
            else:
                break
        
        # populate attr with values from remaining columns, if present
        for i in range(num_bed_columns,len(items)):
            name, formatter = column_formatters[i-num_bed_columns] 
            attr[name] = formatter(items[i])
        
        # stash order of columns for export
        if num_bed_columns > 0:
            attr["_bedx_column_order"] = [X[0] for X in column_formatters]
    
        # convert color to hex string
        try:
            attr["color"] = get_str_from_rgb255(tuple([int(X) for X in attr["color"].split(",")]))
        except ValueError:
            attr["color"] = "#000000"
    
        # sanity check on thickstart and thickend
        if attr["thickstart"] == attr["thickend"]: # if coding region is 0 length, RNA is non-coding
            attr["thickstart"] = attr["thickend"] = chrom_start
        elif any([attr["thickstart"] is None, attr["thickend"] is None]):
            attr["thickstart"] = attr["thickend"] = chrom_start
        elif attr["thickstart"] < 0 or attr["thickend"] < 0:
            attr["thickstart"] = attr["thickend"] = chrom_start
        
        # convert blocks to GenomicSegments
        num_frags    = int(attr["blocks"])
        frag_sizes   = [int(X) for X in attr["blocksizes"].strip(",").split(",")[:num_frags]]
        frag_offsets = [int(X) for X in attr["blockstarts"].strip(",").split(",")[:num_frags]]
        for i in range(0,num_frags):
            frag_start = chrom_start + frag_offsets[i]
            frag_end   = frag_start  + frag_sizes[i]
            frags.append(GenomicSegment(chrom,frag_start,frag_end,strand))

        # clean up attr
        for k in ("blocks","blocksizes","blockstarts"):
            attr.pop(k)
    
        return SegmentChain(*frags,**attr)
    
    @staticmethod
    def from_psl(psl_line):
        """Create a |SegmentChain| from a line from a `PSL`_ (BLAT) file
        
        Parameters
        ----------
        psl_line : str
            Line from a `PSL`_ file
        
        See the `PSL spec <http://pombe.nci.nih.gov/genome/goldenPath/help/blatSpec.html>`_
        
        Returns
        -------
        |SegmentChain|
        """
        cdef:
            list items = psl_line.strip().split("\t")        
            list segs = []
            list block_starts, q_starts, t_starts
            dict attr = {}
            long t_start, block_size
            GenomicSegment seg

        attr["type"]             = "alignment"
        attr["query_name"]       = items[9]
        attr["match_length"]     = int(items[0])
        attr["mismatches"]       = int(items[1])
        attr["rep_matches"]      = int(items[2])
        attr["N"]                = int(items[3])
        attr["query_gap_count"]  = int(items[4])
        attr["query_gap_bases"]  = int(items[5])
        attr["target_gap_count"] = int(items[6])
        attr["target_gap_bases"] = int(items[7])
        attr["strand"]           = items[8]
        attr["query_length"]     = int(items[10])
        attr["query_start"]      = int(items[11])
        attr["query_end"]        = int(items[12])
        attr["target_name"]      = items[13]
        attr["target_length"]    = int(items[14])
        attr["target_start"]     = int(items[15])
        attr["target_end"]       = int(items[16])
        attr["ID"]               = attr["query_name"]
        #block_count           = int(items[17])

        block_sizes = [int(X) for X in items[18].strip(",").split(",")]
        q_starts    = [int(X) for X in items[19].strip(",").split(",")]
        t_starts    = [int(X) for X in items[20].strip(",").split(",")]        

        attr["q_starts"] = q_starts
        attr["t_starts"] = t_starts
        
        for t_start, block_size in zip(t_starts,block_sizes):
            seg = GenomicSegment(attr["target_name"],
                                 t_start,
                                 t_start + block_size,
                                 attr["strand"])
            segs.append(seg)
        
        return SegmentChain(*segs,**attr)        

    def __copy__(self): 
        return SegmentChain(*self,**self.attr)

    def __deepcopy__(self,memo):
        new_ivs  = copy.deepcopy(self._segments)
        new_attr = copy.deepcopy(self.attr)
        return SegmentChain(*new_ivs,**new_attr)    



# TODO: fix properties
cdef class Transcript(SegmentChain):
    """Subclass of |SegmentChain| specifically for transcripts.
    In addition to coordinate-conversion, count fetching, sequence fetching,
    and various other methods inherited from |SegmentChain|, |Transcript|
    provides convenience methods for fetching sub-chains corresponding to 
    CDS features, 5' UTRs, and 3' UTRs.

    Attributes
    ----------
    cds_genome_start : int or None
        Leftmost position in genomic coordinates of coding region, 0-indexed

    cds_genome_end : int or None
        Rightmost position in genomic coordinates of coding region, 0-indexed
        and half-open

    cds_start : int or None
        Stranded position relative to 5' end of transcript at which coding region starts
        (note: for minus-strand features this will be higher in genomic
        coordinates than `cds_end`).

    cds_end : int or None
        Stranded position relative to 5' end of transcript at which coding region ends
        (note: for minus-strand features this will be lower in genomic coordinates
        than `cds_start`).

    iv : |GenomicSegment|
        A GenomicSegment spanning the endpoints of the Transcript

    strand : str
        The chromosome strand ('+', '-', or '.')

    chrom : str
        The chromosome name

    attr : dict
        Miscellaneous attributes
    """
    
    def __init__(self,*ivs,**attr):
        """Create a |Transcript|
        
        Parameters
        ----------
        *ivs : |GenomicSegment|
            0 or more |GenomicSegments| (exons)

        **attr : dict
            keyword attributes

        attr["cds_genome_start"] : int or None
            genome coordinate of CDS start, if any
                         
        attr["cds_genome_end"] : int or None
            genome coordinate of CDS end, if any
    
        attr["type"] : str
            If provided, a feature type used for GTF2/GFF3 export
            Otherwise, set to "mRNA"
        
        attr["ID"] : str
            If provided, a unique ID for the |Transcript|.
            Otherwise, generated from genomic coordinates
        
        attr["transcript_id"] : str
            If provided, a transcript_id used for `GTF2`_ export.
            Otherwise, generated from genomic coordinates.
        
        attr["gene_id"] : str
            If provided, a gene_id used for `GTF2`_ export
            Otherwise, generated from genomic coordinates.
        """
        SegmentChain.__init__(self,*ivs,**attr)
        if "type" not in attr:
            attr["type"] = "mRNA"
            
        self.cds_genome_start = long(attr.get("cds_genome_start",None))
        self.cds_genome_end   = long(attr.get("cds_genome_end",None))
        self._update()
 
    property cds_start:
        def __get__(self):
            return self.cds_start
        def __set__(self, long val):
            self.cds_start = val
            if val is not None:
                if self.strand == "-":
                    self.cds_genome_end = self._position_hash[val+1]
                else:
                    self.cds_genome_start = self._position_hash[val]
            else:
                 if self.strand == "-":
                    self.cds_genome_end = None
                 else:
                    self.cds_genome_start = None

    property cds_end:
        def __get__(self):
            return self.cds_end
        def __set__(self, long val):
            self.cds_end = val
            if val is not None:
                if self.strand == "-":
                    self.cds_genome_start = self._position_hash[val]
                else:
                    self.cds_genome_end = self._position_hash[val]
            else:
                if self.strand == "-":
                    self.cds_genome_start = None
                else:
                    self.cds_genome_end = None

    property cds_genome_start:
        def __get__(self):
            return self.cds_genome_start
        def __set__(self, long val):
            self.cds_genome_start = val
            if val is not None:
                if self.strand == "-":
                    self.cds_end = self._inverse_hash[val] + 1
                else:
                    self.cds_start = self._inverse_hash[val]
            else:
                if self.strand == "-":
                    self.cds_end = None
                else:
                    self.cds_start = None

    property cds_genome_end:
        def __get__(self):
            return self.cds_genome_start
        def __set__(self, long val):
            self.cds_genome_end = val
            if val is not None:
                if self.strand == "-":
                    self.cds_start = self._inverse_hash[val] - 1
                else:
                    self.cds_end = self._inverse_hash[val]
            else:
                if self.strand == "-":
                    self.cds_start = None
                else:
                    self.cds_end = None

    cpdef void _update(self):
        SegmentChain._update(self)
        if self.cds_genome_start is not None:
            if self.strand == "+":
                self.cds_start = self.get_segmentchain_coordinate(self.chrom,self.cds_genome_start,self.strand,stranded=True)
                
                # this is in a try-catch because if the half-open cds_end coincides
                # with the end of an exon, it will not be in the end-inclusive position
                try:
                    self.cds_end = self.get_segmentchain_coordinate(self.chrom,self.cds_genome_end, self.strand,stranded=True)
                except KeyError:
                    # minus one, plus one corrections because end-exclusive genome
                    # position will not be in position hash if it coincides with
                    # the end of any exon
                    self.cds_end   = 1 + self.get_segmentchain_coordinate(self.chrom,self.cds_genome_end - 1, self.strand,stranded=True)
            else:
                # likewise for minus-strand
                # both this adjustment and the one above for plus-strand features
                # have been thoroughly tested by examining BED files exported
                # for this purpose
                self.cds_start = self.get_segmentchain_coordinate(self.chrom,self.cds_genome_end - 1, self.strand,stranded=True)
                self.cds_end   = 1 + self.get_segmentchain_coordinate(self.chrom,self.cds_genome_start,self.strand,stranded=True)
        else:
            self.cds_start = None
            self.cds_end   = None


    def get_name(self):
        """Return the name of `self`, first searching through
        `self.attr` for the keys `transcript_id`, `ID`, `Name`, and `name`.
        If no value is found, :meth:`Transcript.__str__` is used.
        
        Returns
        -------
        str
            Returns in order of preference, `transcript_id`, `ID`, `Name`,
            or `name` from `self.attr`. If not found, returns ``str(self)``
        """
        name = self.attr.get("transcript_id",
               self.attr.get("ID",
               self.attr.get("Name",
               self.attr.get("name",
                             str(self)))))
        return name
    
    def get_cds(self,**extra_attr):
        """Retrieve |SegmentChain| covering the coding region of `self`, including the stop codon.
        If no coding region is present, returns an empty |SegmentChain|.
        
        The following attributes are passed from `self.attr` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`
            #. ID, generated as `"%s_CDS % self.get_name()`


        Parameters
        ----------
        extra_attr : keyword arguments
            Values that will be included in the CDS subchain's `attr` dict.
            These can be used to overwrite values already present.
        
        Returns
        -------
        |SegmentChain|
            CDS region of `self` if present, otherwise empty |SegmentChain|
        """
        cdef SegmentChain my_segmentchain
        if self.cds_genome_start is not None and self.cds_genome_end is not None:
            my_segmentchain = self.get_subchain(self.cds_start,
                                                self.cds_end,
                                                stranded=True,
                                                thickstart=self.cds_genome_start,
                                                thickend=self.cds_genome_end,
                                                cds_start=self.cds_start,
                                                cds_end=self.cds_end,
                                                cds_genome_start=self.cds_genome_start,
                                                cds_genome_end=self.cds_genome_end,
                                                ID="%s_CDS" % self.get_name(),
                                                transcript_id=self.get_name(),
                                                type="CDS",gene_id=self.get_gene())
            my_segmentchain.attr.update(extra_attr)
            
        return my_segmentchain   
    
    def get_utr5(self,**extra_attr):
        """Retrieve sub-|SegmentChain| covering 5'UTR of `self`.
        If no coding region, returns an empty |SegmentChain|

        The following attributes are passed from `self.attr` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`
            #. ID, generated as `"%s_5UTR" % self.get_name()`


        Parameters
        ----------
        extra_attr : keyword arguments
            Values that will be included in the 5'UTR subchain's `attr` dict.
            These can be used to overwrite values already present.


        Returns
        -------
        |SegmentChain|
            5' UTR region of `self` if present, otherwise empty |SegmentChain|
        """
        cdef:
            SegmentChain my_segmentchain = SegmentChain()
            dict attr = {}
        if self.cds_genome_start is not None and self.cds_genome_end is not None:
            my_segmentchain = self.get_subchain(0,self.cds_start,stranded=True)
            attr["type"] = "5UTR"
            attr["gene_id"] = self.get_gene()
            attr["transcript_id"] = self.get_gene()
            attr["ID"] = "%s_5UTR" % self.get_name()
            my_segmentchain.attr = attr
            my_segmentchain.attr.update(extra_attr)

        return my_segmentchain   
    
    def get_utr3(self,**extra_attr):
        """Retrieve sub-|SegmentChain| covering 3'UTR of `self`, excluding
        the stop codon. If no coding region, returns an empty |SegmentChain|
        
        The following attributes are passed from ``self.attr`` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`
            #. ID, generated as `"%s_3UTR" % self.get_name()`


        Parameters
        ----------
        extra_attr : keyword arguments
            Values that will be included in the 3' UTR subchain's `attr` dict.
            These can be used to overwrite values already present.


        Returns
        -------
        |SegmentChain|
            3' UTR region of `self` if present, otherwise empty |SegmentChain|
        """
        cdef:
            SegmentChain my_segmentchain
            dict attr = {}

        if self.cds_genome_start is not None and self.cds_genome_end is not None:
            my_segmentchain = self.get_subchain(self.cds_end,self.length,stranded=True)
            attr["type"] = "3UTR"
            attr["gene_id"] = self.get_gene()
            attr["transcript_id"] = self.get_gene()
            attr["ID"] = "%s_3UTR" % self.get_name()
            my_segmentchain.attr.update(attr)
            my_segmentchain.attr.update(extra_attr)
            
        return my_segmentchain   

    def as_gtf(self, str feature_type="exon", bint escape=True, list excludes=[]):
        """Format `self` as a `GTF2`_ block. |GenomicSegments| are formatted
        as `GTF2`_ `'exon'` features. Coding regions, if peresent, are formatted
        as `GTF2`_ `'CDS'` features. Stop codons are excluded in the `'CDS'` features,
        per the `GTF2`_ specification, and exported separately.

        All attributes from `self.attr` are propagated to the exon and CDS
        features that are generated.

         
        Parameters
        ----------
        feature_type : str
            If not None, overrides the `'type'` attribute of `self.attr`
        
        escape : bool, optional
            URL escape tokens in column 9 of `GTF`_ output (Default: `True`)
        
        
        Returns
        -------
        str
            Block of GTF2-formatted text


        Notes
        -----
        `gene_id` and `transcript_id` are required
            The `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_ requires
            that attributes `gene_id` and `transcript_id` be defined. If these
            are not present in `self.attr`, their values will be guessed 
            following the rules in :py:meth:`SegmentChain.get_gene` and 
            :py:meth:`SegmentChain.get_name`, respectively.
        
        Beware of attribute loss
            To save memory, only the attributes shared by all of the individual
            sub-features (e.g. exons) that were used to assemble this |Transcript|
            have been stored in `self.attr`. This means that upon re-export to `GTF2`_,
            these sub-features will be lacking any attributes that were specific
            to them individually. Formally, this is compliant with the 
            `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_, which states
            explicitly that only the attributes `gene_id` and `transcript_id`
            are supported.
            
        Columns of `GTF2`_ are as follows
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if coding) 
                9     Attributes. "gene_id" and "transcript_id" are required                        
            ======== =========
        
        For more info
            - `GTF2 file format specification <http://mblab.wustl.edu/GTF22.html>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_        
        """
        cdef:
            str stmp
            SegmentChain cds_chain_temp, start_codon_chain, stop_codon_chain
            dict child_chain_attr
            numpy.ndarray[LONG_t,ndim=1] cds_positions
            GenomicSegment span = self.spanning_segment
            my_strand = span.strand
            my_chrom = span.chrom

        stmp  = SegmentChain.as_gtf(self,feature_type=feature_type,escape=escape,excludes=[])
        cds_chain_temp = self.get_cds()
        if len(cds_chain_temp) > 0:
            child_chain_attr  = copy.deepcopy(self.attr)
            child_chain_attr.pop("type")
            cds_positions = cds_chain_temp.get_position_list()
            
            # remove stop codons from CDS, per GTF2 spec
            if my_strand == "+":
                cds_positions = cds_positions[:-3]
            else:
                cds_positions = cds_positions[3:]
            cds_chain = SegmentChain(*positionlist_to_segments(my_chrom,
                                                               my_strand,
                                                               cds_positions),
                                   type="CDS",**child_chain_attr)
            stmp += cds_chain.as_gtf(feature_type="CDS",escape=escape,excludes=excludes)
            
            start_codon_chain = cds_chain_temp.get_subchain(0, 3)
            start_codon_chain.attr.update(child_chain_attr)
            stmp += start_codon_chain.as_gtf(feature_type="start_codon",escape=escape,excludes=excludes)
    
            stop_codon_chain = cds_chain_temp.get_subchain(cds_chain_temp.get_length()-3, cds_chain_temp.get_length())
            stop_codon_chain.attr.update(child_chain_attr)
            stmp += stop_codon_chain.as_gtf(feature_type="stop_codon",escape=escape,excludes=excludes)

        return stmp

    def as_gff3(self, bint escape=True, list excludes=[], str rna_type="mRNA"):
        """Format a |Transcript| as a block of `GFF3`_ output, following
        the schema set out in the `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_
        
        The |Transcript| will be formatted according to the following rules:
        
          1. A feature of type `rna_type` will be created, with `Parent` attribute
             set to the value of ``self.get_gene()``, and `ID` attribute
             set to ``self.get_name()``
        
          2. For each |GenomicSegment| in `self`, a child feature of type
             `exon` will be created. The `Parent` attribute of these features
             will be set to the value of ``self.get_name()``. These will
             have unique IDs generated from ``self.get_name()``.

          3. If `self` is coding (i.e. has none-`None` value for
             `self.cds_genome_start` and `self.cds_genome_end`), child features
             of type `'five_prime_UTR'`, `'CDS'`, and `'three_prime_UTR'` will be created,
             with `Parent` attributes set to ``self.get_name()``. These will
             have unique IDs generated from ``self.get_name()``.
        
        
        Parameters
        ----------
        escape : bool, optional
            Escape tokens in column 9 of `GFF3`_ output (Default: `True`)
        
        excludes : list, optional
            List of attribute key names to exclude from column 9
            (Default: `[]`)
        
        rna_type : str, optional
            Feature type to export RNA as (e.g. `'tRNA'`, `'noncoding_RNA'`,
            et c. Default: `'mRNA'`)

        
        Returns
        -------
        str
            Multiline block of `GFF3`_-formatted text


        Notes
        -----
        Beware of attribute loss
            This |Transcript| was assembled from multiple individual component
            features (e.g. single exons), which may or may not have had their own
            unique attributes in their original annotation. To reduce overhead, 
            these individual attributes (if they were present) have not been
            (entirely) stored, and consequently will not (all) be exported.
            If this poses problems, consider instead importing, modifying, and
            exporting the component features

        GFF3 schemas vary
            Different GFF3s have different schemas (parent-child relationships
            between features). Here we adopt the commonly-used schema set by
            `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_,
            which may or may not match your schema.

        Columns of `GFF3`_ are as follows
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if coding) 
                9     Attributes                       
            ======== =========

        For futher information, see
            - `GFF3 file format specification <http://www.sequenceontology.org/gff3.shtml>`_
            - `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`
            - `SO releases <http://sourceforge.net/projects/song/files/SO_Feature_Annotation/>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        """
        cdef:
            str gene_id       = self.get_gene()
            str transcript_id = self.get_name()
            str ftype, chain, my_id
            list ltmp = []
            list parts = []
            dict child_attr = copy.deepcopy(self.attr)
            SegmentChain feature
            GenomicSegment iv
            int n

#        keys_to_pop = ("ID",)
#        for k in keys_to_pop:
#            if k in child_attr:
#                child_attr.pop(k)

        # mRNA feature
        feature = SegmentChain(self.spanning_segment,ID=transcript_id,Parent=gene_id,type=rna_type)
        ltmp.append(feature.as_gff3(excludes=excludes,escape=escape))

        # child features
        child_attr["Parent"] = transcript_id

        # exon feature
        child_attr["type"] = "exon"
        for n, seg in enumerate(self):
            my_id   = "%s:exon:%s" % (transcript_id,n)
            child_attr["ID"]   = my_id
            feature = SegmentChain(seg,**child_attr)
            ltmp.append(feature.as_gff3(excludes=excludes,escape=escape))
        
        # CDS & UTRs
        if self.cds_genome_start is not None:
            parts = [("five_prime_UTR", self.get_utr5()),
                     ("CDS",            self.get_cds()),
                     ("three_prime_UTR",self.get_utr3()),
                    ]
            for ftype, chain in parts:
                child_attr["type"]   = ftype
                child_attr["Parent"] = transcript_id
                for n, seg in enumerate(self):
                    my_id   = "%s:%s:%s" % (transcript_id,ftype,n)
                    child_attr["ID"]   = my_id
                    feature = SegmentChain(seg,**child_attr)
                    ltmp.append(feature.as_gff3(excludes=excludes,escape=escape))

        return "".join(ltmp)
    
    def as_bed(self,as_int=True,color=None,extra_columns=None):
        """Format `self` as a BED12[+X] line, assigning CDS boundaries 
        to the thickstart and thickend columns from `self.attr`

        If the |SegmentChain| was imported as a `BED`_ file with extra columns,
        these will be output in the same order, after the `BED`_ columns.
        
        Parameters
        ----------
        as_int : bool, optional
            Force "score" to integer (Default: True)
    
        color : str or None, optional
            Color represented as RGB hex string.
            If not none, overrides the color in `self.attr["color"]`
    
        extra_columns : None or list, optional
            If `None`, and the |Transcript| was imported using the `extra_columns`
            keyword of :meth:`~yeti.genomics.roitools.Transcript.from_bed`,
            the |Transcript| will be exported in BED 12+X format, in which
            extra columns are in the same order as they were upon import. If no extra columns
            were present, the |Transcript| will be exported a aa BED12 line.

            If a list of attribute names, these attributes will be exported as
            extra columns in order, overriding whatever happened upon import. 
            If an attribute name is not in the `attr` dict of the |Transcript|,
            it will be exported with the default empty value "".

            If an empty list, no extra columns will be exported; the |Transcript|
            will be formatted as a BED12 line.

    
        Returns
        -------
        str
            Line of BED12-formatted text
            
        
        Notes
        -----
        BED12 columns are as follows
            ======== =========
            Column   Contains
            ======== =========
               0     Contig or chromosome
               1     Start of first block in feature (0-indexed)
               2     End of last block in feature (half-open)
               3     Feature name
               4     Feature score
               5     Strand
               6     thickstart
               7     thickend
               8     Feature color as RGB tuple
               9     Number of blocks in feature
               10    Block lengths
               11    Block starts, relative to start of first block
            ======== =========

        Fore more information
            See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        """
        return SegmentChain.as_bed(self,
                                   thickstart=self.cds_genome_start,
                                   thickend=self.cds_genome_end,
                                   as_int=as_int,
                                   color=color,
                                   extra_columns=extra_columns)

    @staticmethod
    def from_bed(str line,extra_columns=0):
        """Create a |Transcript| from a BED line with 4 or more columns.
        `thickstart` and `thickend` columns, if present, are assumed to specify
        CDS boundaries, a convention that, while common, is formally outside the
        `BED`_ specification.
    
    	See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
    	for more details.

        Parameters
        ----------
        line
            Line from a BED file with at least 4 columns

        extra_columns: int or list, optional
            Extra, non-BED columns in `BED`_ file corresponding to feature
            attributes. This is common in `ENCODE`_-specific `BED`_ variants.
            
            if `extra-columns` is:
            
              - an :class:`int`: it is taken to be the
                number of attribute columns. Attributes will be stored in
                the `attr` dictionary of the |SegmentChain|, under names like
                `custom0`, `custom1`, ... , `customN`.

              - a :class:`list` of :class:`str`, it is taken to be the names
                of the attribute columns, in order, from left to right in the file.
                In this case, attributes in extra columns will be stored under

              - a :class:`list` of :class:`tuple`, each tuple is taken
                to be a pair of `(attribute_name, formatter_func)`. In this case,
                the value of `attribute_name` in the `attr` dict of the |SegmentChain|
                will be set to `formatter_func(column_value)`.
            
            (Default: 0)
                
    
        Returns
        -------
        |Transcript|
        """
        cdef:
            SegmentChain segchain = SegmentChain.from_bed(line,extra_columns=extra_columns)
            segments = segchain._segments
            attr = segchain.attr
            Transcript transcript

        attr.pop("type") # default type for SegmentChain is "exon". We want to use "mRNA"
        attr["cds_genome_start"] = attr["thickstart"]
        attr["cds_genome_end"]   = attr["thickend"]
        attr.pop("thickstart")
        attr.pop("thickend")
        if attr["cds_genome_start"] == attr["cds_genome_end"]:
            attr["cds_genome_start"] = attr["cds_genome_end"] = None

        transcript = Transcript(*segments,**attr)
    
        return transcript
    
    @staticmethod
    def from_psl(str psl_line):
        cdef:
            SegmentChain segchain = SegmentChain.from_psl(psl_line)
            dict attr = segchain.attr

        attr["cds_genome_start"] = None
        attr["cds_genome_end"] = None
        return Transcript(*segchain,**attr)
