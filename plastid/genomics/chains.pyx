"""This module contains experimental Cython implementations of data structures 
in plastid. The implementation of SegmenChain below works and will be adopted
after further tets. The implementation of Transcript below does not, and
should be avoided for now
"""
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

from plastid.util.services.exceptions import DataWarning
from plastid.util.services.colors import get_str_from_rgb255, get_rgb255_from_str
from plastid.readers.gff_tokens import make_GFF3_tokens, \
                                    make_GTF2_tokens


from plastid.genomics.c_roitools cimport GenomicSegment, \
                                      positionlist_to_segments, \
                                      positions_to_segments, \
                                      strand_to_str, str_to_strand
from plastid.genomics.c_roitools import NullSegment
from plastid.genomics.c_common cimport ExBool, true, false, bool_exception, \
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
# higher-order classes that handle multi-feature structures, like transcripts
# or alignments
#===============================================================================

cdef class s1(object):

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
        cdef:
            int num_segs = len(self)
            list segs
            GenomicSegment seg0
        
        self.sort()
        self._get_position_hash()

        self.reset_masks()
       
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
        self.length = length

    cpdef void sort(self):
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
    
    def get_position_list(self,copy = False):
        return self.c_get_position_list(copy)

    cdef numpy.ndarray c_get_position_list(self, bint copy):
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._position_hash,dtype=long) 
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cpdef set get_position_set(self):
        return set(self.get_position_list())
    
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
       
    cdef tuple check_segments(self, tuple segments):
        cdef:
            GenomicSegment seg, seg0
            GenomicSegment span = self.spanning_segment
            str my_chrom  = span.chrom
            Strand my_strand = span.c_strand
            int i = 0
            int length = len(segments)
            bint strandprob = False
            bint chromprob = False

        if length > 0:
            seg0 = segments[0]
            if len(self._segments) == 0:
                my_chrom  = seg0.chrom
                my_strand = seg0.c_strand

            while i < length:
                seg = segments[i]
                if seg.chrom != my_chrom:
                    chromprob = True
                    break
                    
                if seg.c_strand != my_strand:
                    strandprob = True
                    break

                i += 1

            if strandprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))
            if chromprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
        return (my_chrom,strand_to_str(my_strand))

    cdef c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            list positions = list(self.c_get_position_list(False))
            GenomicSegment seg
            int length = len(segments)

        if length > 0:
            my_chrom, my_strand = self.check_segments(segments)
            # add new positions
            for seg in segments:
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positionlist_to_segments(my_chrom,my_strand,sorted(list(set(positions))))
            self._update()

    def add_segments(self,*segments):
        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
    
    cpdef void reset_masks(self):
        cdef:
            long length = self.length
            int [:] pmask = array.clone(mask_template,length,True)

        self._position_mask = pmask
        self._mask_segments = []
        self.masked_length = length


cdef tuple check_segments(s2 chain, tuple segments):
    cdef:
        GenomicSegment seg, seg0
        GenomicSegment span = chain.spanning_segment
        str my_chrom  = span.chrom
        Strand my_strand = span.c_strand
        int i = 0
        int length = len(segments)
        bint strandprob = False
        bint chromprob = False

    if length > 0:
        seg0 = segments[0]
        if len(chain._segments) == 0:
            my_chrom  = seg0.chrom
            my_strand = seg0.c_strand

        while i < length:
            seg = segments[i]
            if seg.chrom != my_chrom:
                chromprob = True
                break
                
            if seg.c_strand != my_strand:
                strandprob = True
                break

            i += 1

        if strandprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(chain)))
        if chromprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(chain)))
    return (my_chrom,strand_to_str(my_strand))

cdef void add_segments(s2 chain, tuple segments):
    cdef:
        str my_chrom, my_strand
        GenomicSegment seg, seg0
        list positions = []
        int length = len(segments)
        int num_segs
        list new_segments

        cdef long total_length
        cdef long [:] my_hash
        cdef long x, c
        cdef dict ihash = {}

    if length > 0:
        my_chrom, my_strand = check_segments(chain,segments)
        # add new positions
        for seg in segments:
            positions.extend(range(seg.start,seg.end))
        
        # reset variables
        new_segments = positionlist_to_segments(my_chrom,my_strand,sorted(list(set(positions))))
        chain._segments = new_segments
        num_segs = len(new_segments)
        
        if num_segs == 0:
            chain.spanning_segment = NullSegment 
        elif num_segs == 1:
            chain.spanning_segment = chain._segments[0]
        elif num_segs >1:
            segs = chain._segments
            seg0 = segs[0]
            chain.spanning_segment = GenomicSegment(seg0.chrom,
                                                   seg0.start,
                                                   segs[-1].end,
                                                   seg0.strand)

        length = sum([len(X) for X in segments])
        my_hash = array.clone(hash_template,length,False)  

        c = 0
        for seg in new_segments:
            x = seg.start
            while x < seg.end:
                my_hash[c] = x
                ihash[x] = c
                c += 1
                x += 1
        
        chain._position_hash = my_hash
        chain._inverse_hash  = ihash
        chain.length = length


# as it is, this is far faster for
#   - init empty (5x)
#   - init with segments (4x)
#
# slower for:
#   - add segments (10%)
#   - reset masks (10x)
cdef class s2(object):

    def __cinit__(self, *segments, **attr):
        self._segments        = []
        self._mask_segments   = []
        self._position_hash   = array.clone(hash_template,0,False)
        self._position_mask   = array.clone(mask_template,0,False)
        self._inverse_hash    = {}
        self.length = 0
        self.masked_length = 0
        self.spanning_segment = NullSegment
        self.attr = attr
        if "type" not in attr:
            attr["type"] = "exon"
        self.attr = attr
        add_segments(self,segments)

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
        cdef:
            int num_segs = len(self)
            list segs
            GenomicSegment seg0
        
        self.sort()
        self._get_position_hash()

        self.reset_masks()
       
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
        self.length = length

    cpdef void sort(self):
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
    
    def get_position_list(self,copy = False):
        return self.c_get_position_list(copy)

    cdef numpy.ndarray c_get_position_list(self, bint copy):
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._position_hash,dtype=long) 
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cpdef set get_position_set(self):
        return set(self.get_position_list())
    
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
       
    cdef tuple check_segments(self, tuple segments):
        cdef:
            GenomicSegment seg, seg0
            GenomicSegment span = self.spanning_segment
            str my_chrom  = span.chrom
            Strand my_strand = span.c_strand
            int i = 0
            int length = len(segments)
            bint strandprob = False
            bint chromprob = False

        if length > 0:
            seg0 = segments[0]
            if len(self._segments) == 0:
                my_chrom  = seg0.chrom
                my_strand = seg0.c_strand

            while i < length:
                seg = segments[i]
                if seg.chrom != my_chrom:
                    chromprob = True
                    break
                    
                if seg.c_strand != my_strand:
                    strandprob = True
                    break

                i += 1

            if strandprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))
            if chromprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
        return (my_chrom,strand_to_str(my_strand))

    cdef c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            list positions = list(self.c_get_position_list(False))
            GenomicSegment seg
            int length = len(segments)

        if length > 0:
            my_chrom, my_strand = self.check_segments(segments)
            # add new positions
            for seg in segments:
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positionlist_to_segments(my_chrom,my_strand,sorted(list(set(positions))))
            self._update()

    def add_segments(self,*segments):
        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
    
    cpdef void reset_masks(self):
        cdef:
            long length = self.length
            int [:] pmask = array.clone(mask_template,length,True)

        self._position_mask = pmask
        self._mask_segments = []
        self.masked_length = length


cdef tuple s3_check_segments(s3 chain, tuple segments):
    cdef:
        GenomicSegment seg, seg0
        GenomicSegment span = chain.spanning_segment
        str my_chrom  = span.chrom
        Strand my_strand = span.c_strand
        int i = 0
        int length = len(segments)
        bint strandprob = False
        bint chromprob = False

    if length > 0:
        seg0 = segments[0]
        if len(chain._segments) == 0:
            my_chrom  = seg0.chrom
            my_strand = seg0.c_strand

        while i < length:
            seg = segments[i]
            if seg.chrom != my_chrom:
                chromprob = True
                break
                
            if seg.c_strand != my_strand:
                strandprob = True
                break

            i += 1

        if strandprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(chain)))
        if chromprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(chain)))
    return (my_chrom,strand_to_str(my_strand))


cdef void s3_setup(s3 chain, tuple segments):
    cdef:
        str my_chrom, my_strand
        GenomicSegment seg, seg0
        list positions = []
        int length = len(segments)
        int num_segs
        list new_segments

        cdef long total_length
        cdef long [:] my_hash
        cdef long x, c
        cdef dict ihash = {}
        int [:] pmask

    if length > 0:
        my_chrom, my_strand = s3_check_segments(chain,segments)
        # add new positions
        for seg in segments:
            positions.extend(range(seg.start,seg.end))
        
        # reset variables
        new_segments = positionlist_to_segments(my_chrom,my_strand,sorted(list(set(positions))))
        chain._segments = new_segments
        num_segs = len(new_segments)
        
        if num_segs == 0:
            chain.spanning_segment = NullSegment 
        elif num_segs == 1:
            chain.spanning_segment = new_segments[0]
        elif num_segs >1:
            segs = new_segments
            seg0 = segs[0]
            chain.spanning_segment = GenomicSegment(seg0.chrom,
                                                   seg0.start,
                                                   segs[-1].end,
                                                   seg0.strand)

        length = sum([len(X) for X in segments])
        my_hash = array.clone(hash_template,length,False)  

        c = 0
        for seg in new_segments:
            x = seg.start
            while x < seg.end:
                my_hash[c] = x
                ihash[x] = c
                c += 1
                x += 1
        
        chain._position_hash = my_hash
        chain._inverse_hash  = ihash
        chain.length = length
        pmask = array.clone(mask_template,length,True)
        chain._position_mask = pmask

# faster for
#
# slower for:
cdef class s3(object):

    def __cinit__(self, *segments, **attr):
        self._segments        = []
        self._mask_segments   = []
        self._position_hash   = array.clone(hash_template,0,False)
        self._position_mask   = array.clone(mask_template,0,False)
        self._inverse_hash    = {}
        self.length = 0
        self.masked_length = 0
        self.spanning_segment = NullSegment
        if "type" not in attr:
            attr["type"] = "exon"
        self.attr = attr
        s3_setup(self,segments)

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
        cdef:
            int num_segs = len(self)
            list segs
            GenomicSegment seg0
        
        self.sort()
        self._get_position_hash()

        self.reset_masks()
       
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
        self.length = length

    cpdef void sort(self):
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
    
    def get_position_list(self,copy = False):
        return self.c_get_position_list(copy)

    cdef numpy.ndarray c_get_position_list(self, bint copy):
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._position_hash,dtype=long) 
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cpdef set get_position_set(self):
        return set(self.get_position_list())
    
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
       
    cdef tuple check_segments(self, tuple segments):
        cdef:
            GenomicSegment seg, seg0
            GenomicSegment span = self.spanning_segment
            str my_chrom  = span.chrom
            Strand my_strand = span.c_strand
            int i = 0
            int length = len(segments)
            bint strandprob = False
            bint chromprob = False

        if length > 0:
            seg0 = segments[0]
            if len(self._segments) == 0:
                my_chrom  = seg0.chrom
                my_strand = seg0.c_strand

            while i < length:
                seg = segments[i]
                if seg.chrom != my_chrom:
                    chromprob = True
                    break
                    
                if seg.c_strand != my_strand:
                    strandprob = True
                    break

                i += 1

            if strandprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))
            if chromprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
        return (my_chrom,strand_to_str(my_strand))

    cdef c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            list positions = list(self.c_get_position_list(False))
            GenomicSegment seg
            int length = len(segments)

        if length > 0:
            my_chrom, my_strand = self.check_segments(segments)
            # add new positions
            for seg in segments:
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positionlist_to_segments(my_chrom,my_strand,sorted(list(set(positions))))
            self._update()

    def add_segments(self,*segments):
        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
    
    cpdef void reset_masks(self):
        self._position_mask [:] = 0
        self._mask_segments = []
        self.masked_length = self.length


cdef tuple s4_check_segments(s4 chain, tuple segments):
    cdef:
        GenomicSegment seg, seg0
        GenomicSegment span = chain.spanning_segment
        str my_chrom  = span.chrom
        Strand my_strand = span.c_strand
        int i = 0
        int length = len(segments)
        bint strandprob = False
        bint chromprob = False

    if length > 0:
        seg0 = segments[0]
        if len(chain._segments) == 0:
            my_chrom  = seg0.chrom
            my_strand = seg0.c_strand

        while i < length:
            seg = segments[i]
            if seg.chrom != my_chrom:
                chromprob = True
                break
                
            if seg.c_strand != my_strand:
                strandprob = True
                break

            i += 1

        if strandprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(chain)))
        if chromprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(chain)))
    return (my_chrom,strand_to_str(my_strand))

cdef void s4_setup(s4 chain, tuple segments):
    cdef:
        str my_chrom, my_strand
        GenomicSegment seg, seg0
        list positions = []
        int length = len(segments)
        int num_segs
        list new_segments

        cdef long total_length
        cdef long [:] my_hash
        cdef long x, c
        cdef dict ihash = {}
        int [:] pmask

    if length == 0:
        chain._segments        = []
        chain._mask_segments   = []
        chain._position_hash   = array.clone(hash_template,0,False)
        chain._position_mask   = array.clone(mask_template,0,False)
        chain._inverse_hash    = {}
        chain.length = 0
        chain.masked_length = 0
        chain.spanning_segment = NullSegment
    else:
        my_chrom, my_strand = s4_check_segments(chain,segments)
        # add new positions
        for seg in segments:
            positions.extend(range(seg.start,seg.end))
        
        # reset variables
        new_segments = positionlist_to_segments(my_chrom,my_strand,sorted(list(set(positions))))
        chain._segments = new_segments
        num_segs = len(new_segments)
        
        if num_segs == 0:
            chain.spanning_segment = NullSegment 
        elif num_segs == 1:
            chain.spanning_segment = new_segments[0]
        elif num_segs >1:
            segs = new_segments
            seg0 = segs[0]
            chain.spanning_segment = GenomicSegment(seg0.chrom,
                                                   seg0.start,
                                                   segs[-1].end,
                                                   seg0.strand)

        length = sum([len(X) for X in segments])
        my_hash = array.clone(hash_template,length,False)  

        c = 0
        for seg in new_segments:
            x = seg.start
            while x < seg.end:
                my_hash[c] = x
                ihash[x] = c
                c += 1
                x += 1
        
        chain._position_hash = my_hash
        chain._inverse_hash  = ihash
        chain.length = length
        pmask = array.clone(mask_template,length,True)
        chain._position_mask = pmask

# faster for
#
# slower for:
cdef class s4(object):

    def __cinit__(self, *segments, **attr):
        s4_setup(self,segments)
        if "type" not in attr:
            attr["type"] = "exon"
        self.attr = attr

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
        cdef:
            int num_segs = len(self)
            list segs
            GenomicSegment seg0
        
        self._get_position_hash()
        self.reset_masks()
       
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
        self.length = length

    cpdef void sort(self):
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
    
    def get_position_list(self,copy = False):
        return self.c_get_position_list(copy)

    cdef numpy.ndarray c_get_position_list(self, bint copy):
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._position_hash,dtype=long) 
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cpdef set get_position_set(self):
        return set(self.get_position_list())
    
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
       
    cdef tuple check_segments(self, tuple segments):
        cdef:
            GenomicSegment seg, seg0
            GenomicSegment span = self.spanning_segment
            str my_chrom  = span.chrom
            Strand my_strand = span.c_strand
            int i = 0
            int length = len(segments)
            bint strandprob = False
            bint chromprob = False

        if length > 0:
            seg0 = segments[0]
            if len(self._segments) == 0:
                my_chrom  = seg0.chrom
                my_strand = seg0.c_strand

            while i < length:
                seg = segments[i]
                if seg.chrom != my_chrom:
                    chromprob = True
                    break
                    
                if seg.c_strand != my_strand:
                    strandprob = True
                    break

                i += 1

            if strandprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))
            if chromprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
        return (my_chrom,strand_to_str(my_strand))

    cdef void c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            list positions = list(self.c_get_position_list(False))
            GenomicSegment seg
            int length = len(segments)

        if length > 0:
            my_chrom, my_strand = self.check_segments(segments)
            # add new positions
            for seg in segments:
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positionlist_to_segments(my_chrom,my_strand,sorted(set(positions)))
            self._update()

    def add_segments(self,*segments):
        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
    
    cpdef void reset_masks(self):
        self._position_mask [:] = 0
        self._mask_segments = []
        self.masked_length = self.length


cdef class s5(object):

    def __init__(self,*segments,**attr):
        self._segments        = []
        self._mask_segments   = []
        self._position_hash   = array.clone(hash_template,0,False)
        self._position_mask   = array.clone(mask_template,0,False)
        self._inverse_hash    = {}
        self.length = 0
        self.masked_length = 0
        self.spanning_segment = NullSegment
        if "type" not in attr:
            attr["type"] = "exon"

        self.attr = attr
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
        cdef:
            int num_segs = len(self)
            list segs
            GenomicSegment seg0
        
        self.sort()
        self._get_position_hash()

        self.reset_masks()
       
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
        self.length = length

    cpdef void sort(self):
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
    
    def get_position_list(self,copy = False):
        return self.c_get_position_list(copy)

    cdef numpy.ndarray c_get_position_list(self, bint copy):
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._position_hash,dtype=long) 
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cpdef set get_position_set(self):
        return set(self.get_position_list())
    
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
       
    cdef tuple check_segments(self, tuple segments):
        cdef:
            GenomicSegment seg, seg0
            GenomicSegment span = self.spanning_segment
            str my_chrom  = span.chrom
            Strand my_strand = span.c_strand
            int i = 0
            int length = len(segments)
            bint strandprob = False
            bint chromprob = False

        if length > 0:
            seg0 = segments[0]
            if len(self._segments) == 0:
                my_chrom  = seg0.chrom
                my_strand = seg0.c_strand

            while i < length:
                seg = segments[i]
                if seg.chrom != my_chrom:
                    chromprob = True
                    break
                    
                if seg.c_strand != my_strand:
                    strandprob = True
                    break

                i += 1

            if strandprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))
            if chromprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
        return (my_chrom,strand_to_str(my_strand))

    cdef c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            list positions = list(self.c_get_position_list(False))
            GenomicSegment seg
            int length = len(segments)

        if length > 0:
            my_chrom, my_strand = self.check_segments(segments)
            # add new positions
            for seg in segments:
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positionlist_to_segments(my_chrom,my_strand,sorted(list(set(positions))))
            self._update()

    def add_segments(self,*segments):
        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
    
    cpdef void reset_masks(self):
        cdef:
            long length = self.length
            int [:] pmask = array.clone(mask_template,length,True)

        self._position_mask = pmask
        self._mask_segments = []
        self.masked_length = length

cdef tuple s6_check_segments(s6 chain, tuple segments):
    cdef:
        GenomicSegment seg, seg0
        GenomicSegment span = chain.spanning_segment
        str my_chrom  = span.chrom
        Strand my_strand = span.c_strand
        int i = 0
        int length = len(segments)
        bint strandprob = False
        bint chromprob = False

    if length > 0:
        seg0 = segments[0]
        if len(chain._segments) == 0:
            my_chrom  = seg0.chrom
            my_strand = seg0.c_strand

        while i < length:
            seg = segments[i]
            if seg.chrom != my_chrom:
                chromprob = True
                break
                
            if seg.c_strand != my_strand:
                strandprob = True
                break

            i += 1

        if strandprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(chain)))
        if chromprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(chain)))
    return (my_chrom,strand_to_str(my_strand))

# faster for
#   - init blank (5x)
#   - init with segments (2x)
#   - add_segments (1.51x) = 856us

cdef class s6(object):

    def __cinit__(self, *segments, **attr):
        cdef:
            str my_chrom, my_strand
            GenomicSegment seg, seg0
            list positions = []
            int length = len(segments)
            int num_segs
            list new_segments

            long total_length, x, c
            long [:] my_hash
            dict ihash = {}
            int [:] pmask

        if "type" not in attr:
            attr["type"] = "exon"
        self.attr = attr
        self._mask_segments = []
        self._segments      = []
        self._inverse_hash    = {}
        self.length = 0
        self.masked_length = 0
        self.spanning_segment = NullSegment
        if length == 0:
            self._position_hash   = array.clone(hash_template,0,False)
            self._position_mask   = array.clone(mask_template,0,False)
        else:
            my_chrom, my_strand = s6_check_segments(self,segments)
            # add new positions
            for seg in segments:
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            new_segments = positionlist_to_segments(my_chrom,my_strand,sorted(set(positions)))
            self._segments = new_segments
            num_segs = len(new_segments)
            
            if num_segs == 0:
                self.spanning_segment = NullSegment 
            elif num_segs == 1:
                self.spanning_segment = new_segments[0]
            elif num_segs >1:
                segs = new_segments
                seg0 = segs[0]
                self.spanning_segment = GenomicSegment(seg0.chrom,
                                                       seg0.start,
                                                       segs[-1].end,
                                                       seg0.strand)

            length = sum([len(X) for X in segments])
            my_hash = array.clone(hash_template,length,False)  

            c = 0
            for seg in new_segments:
                x = seg.start
                while x < seg.end:
                    my_hash[c] = x
                    ihash[x] = c
                    c += 1
                    x += 1
            
            self._position_hash = my_hash
            self._inverse_hash  = ihash
            self.length = length
            self.masked_length = length
            pmask = array.clone(mask_template,length,True)
            self._position_mask = pmask

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

    cdef void _update(self):
        cdef:
            list segments = self._segments
            int num_segs = len(segments)
            GenomicSegment segment, seg0
            long length = sum([len(X) for X in segments])
            long [:] my_hash = array.clone(hash_template,length,False)  
            long x, c
            dict ihash = {}

        self.length = length
        self._position_hash = my_hash
        self._position_mask = array.clone(mask_template,length,True)
        self._inverse_hash  = ihash
        self.c_reset_masks()

        if num_segs == 0:
            self.spanning_segment = NullSegment
        else:
            c = 0
            for segment in segments:
                x = segment.start
                while x < segment.end:
                    my_hash[c] = x
                    ihash[x] = c
                    c += 1
                    x += 1
            
            if num_segs == 1:
                self.spanning_segment = self._segments[0]
            else:
                segs = self._segments
                seg0 = segs[0]
                self.spanning_segment = GenomicSegment(seg0.chrom,
                                                       seg0.start,
                                                       segs[-1].end,
                                                       seg0.strand)

    cpdef void sort(self):
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
    
    def get_position_list2(self):
        return self._position_hash.tolist()

    def get_position_list(self):
        return self.c_get_position_list(False).tolist()

    cdef numpy.ndarray c_get_position_list(self, bint copy):
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._position_hash,dtype=long) 
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cdef set c_get_position_set(self):
        return set(self.c_get_position_list(False).tolist())
    
    def get_position_set(self):
        return self.c_get_position_set()

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
       
    cdef tuple check_segments(self, tuple segments):
        cdef:
            GenomicSegment seg, seg0
            GenomicSegment span = self.spanning_segment
            str my_chrom  = span.chrom
            Strand my_strand = span.c_strand
            int i = 0
            int length = len(segments)
            bint strandprob = False
            bint chromprob = False

        if length > 0:
            seg0 = segments[0]
            if len(self._segments) == 0:
                my_chrom  = seg0.chrom
                my_strand = seg0.c_strand

            while i < length:
                seg = segments[i]
                if seg.chrom != my_chrom:
                    chromprob = True
                    break
                    
                if seg.c_strand != my_strand:
                    strandprob = True
                    break

                i += 1

            if strandprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))
            if chromprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
        return (my_chrom,strand_to_str(my_strand))

    cdef void c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            list positions = self.c_get_position_list(False).tolist()
            GenomicSegment seg
            int length = len(segments)

        if length > 0:
            my_chrom, my_strand = s6_check_segments(self,segments)
            # add new positions
            for seg in segments: # changed
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positionlist_to_segments(my_chrom,my_strand,sorted(set(positions)))
            self._update()

    def add_segments(self,*segments):
        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
    
    cdef void c_reset_masks(self):
        self._position_mask [:] = 0
        self._mask_segments = []
        self.masked_length = self.length

    def reset_masks(self):
        self.c_reset_masks()

cdef tuple s7_check_segments(s7 chain, tuple segments):
    cdef:
        GenomicSegment seg, seg0
        GenomicSegment span = chain.spanning_segment
        str my_chrom  = span.chrom
        Strand my_strand = span.c_strand
        int i = 0
        int length = len(segments)
        bint strandprob = False
        bint chromprob = False

    if length > 0:
        seg0 = segments[0]
        if len(chain._segments) == 0:
            my_chrom  = seg0.chrom
            my_strand = seg0.c_strand

        while i < length:
            seg = segments[i]
            if seg.chrom != my_chrom:
                chromprob = True
                break
                
            if seg.c_strand != my_strand:
                strandprob = True
                break

            i += 1

        if strandprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(chain)))
        if chromprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(chain)))
    return (my_chrom,strand_to_str(my_strand))


cdef class s7(object):

    def __cinit__(self, *segments, **attr):
        cdef:
            str my_chrom, my_strand
            GenomicSegment seg, seg0
            list positions = []
            int length = len(segments)
            int num_segs
            list new_segments
            long total_length
            long c,x
            dict ihash = {}
            long [:] my_hash
            int [:] pmask

        if "type" not in attr:
            attr["type"] = "exon"
        self.attr = attr
        self._mask_segments = []
        self._segments      = []
        self._inverse_hash    = {}
        self.length = 0
        self.masked_length = 0
        self.spanning_segment = NullSegment
        if length == 0:
            self._position_hash   = array.clone(hash_template,0,False)
            self._position_mask   = array.clone(mask_template,0,False)
        else:
            my_chrom, my_strand = s7_check_segments(self,segments)

            # add new positions
            for seg in segments:
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            new_segments = positionlist_to_segments(my_chrom,my_strand,sorted(set(positions)))
            self._segments = new_segments
            num_segs = len(new_segments)
            
            if num_segs == 0:
                self.spanning_segment = NullSegment 
            elif num_segs == 1:
                self.spanning_segment = new_segments[0]
            elif num_segs >1:
                segs = new_segments
                seg0 = segs[0]
                self.spanning_segment = GenomicSegment(seg0.chrom,
                                                       seg0.start,
                                                       segs[-1].end,
                                                       seg0.strand)

            length = sum([len(X) for X in segments])
            my_hash = array.clone(hash_template,length,False)  

            c = 0
            for seg in new_segments:
                x = seg.start
                while x < seg.end:
                    my_hash[c] = x
                    ihash[x] = c
                    c += 1
                    x += 1
            
            self._position_hash = my_hash
            self._inverse_hash  = ihash
            self.length = length
            self.masked_length = length
            pmask = array.clone(mask_template,length,True)
            self._position_mask = pmask

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

    cdef void _update(self):
        cdef:
            list segments = self._segments
            int num_segs = len(segments)
            GenomicSegment segment, seg0
            long length = sum([len(X) for X in segments])
            long x, c
            array.array my_hash = array.clone(hash_template,length,False)  
            dict ihash = {}

        self.length = length
        self._position_hash = my_hash
        self._position_mask = array.clone(mask_template,length,True)
        self._inverse_hash  = ihash
        self.c_reset_masks()

        if num_segs == 0:
            self.spanning_segment = NullSegment
        else:
            c = 0
            for segment in segments:
                x = segment.start
                while x < segment.end:
                    my_hash[c] = x
                    ihash[x] = c
                    c += 1
                    x += 1
            
            if num_segs == 1:
                self.spanning_segment = self._segments[0]
            else:
                segs = self._segments
                seg0 = segs[0]
                self.spanning_segment = GenomicSegment(seg0.chrom,
                                                       seg0.start,
                                                       segs[-1].end,
                                                       seg0.strand)

    cpdef void sort(self):
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
    
    def get_position_list2(self):
        return self._position_hash.aslongs()

    def get_position_list(self):
        return self._position_hash.base.tolist()

    cdef list c_get_position_list(self, bint copy):
        cdef list positions = self._position_hash.base.tolist()
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cdef set c_get_position_set(self):
        return set(self.c_get_position_list(False).tolist())
    
    def get_position_set(self):
        return self.c_get_position_set()

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
       
    cdef tuple check_segments(self, tuple segments):
        cdef:
            GenomicSegment seg, seg0
            GenomicSegment span = self.spanning_segment
            str my_chrom  = span.chrom
            Strand my_strand = span.c_strand
            int i = 0
            int length = len(segments)
            bint strandprob = False
            bint chromprob = False

        if length > 0:
            seg0 = segments[0]
            if len(self._segments) == 0:
                my_chrom  = seg0.chrom
                my_strand = seg0.c_strand

            while i < length:
                seg = segments[i]
                if seg.chrom != my_chrom:
                    chromprob = True
                    break
                    
                if seg.c_strand != my_strand:
                    strandprob = True
                    break

                i += 1

            if strandprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))
            if chromprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
        return (my_chrom,strand_to_str(my_strand))

    cdef void c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            list positions = self.c_get_position_list(False)
            GenomicSegment seg
            int length = len(segments)

        if length > 0:
            my_chrom, my_strand = s7_check_segments(self,segments)
            # add new positions
            for seg in segments: # changed
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positionlist_to_segments(my_chrom,my_strand,sorted(set(positions)))
            self._update()

    def add_segments(self,*segments):
        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
    
    cdef void c_reset_masks(self):
        cdef int [:] pmask = self._position_mask
        pmask [:] = 0
        self._mask_segments = []
        self.masked_length = self.length

    def reset_masks(self):
        self.c_reset_masks()






cdef tuple s8_check_segments(s8 chain, tuple segments):
    cdef:
        GenomicSegment seg, seg0
        GenomicSegment span = chain.spanning_segment
        str my_chrom  = span.chrom
        Strand my_strand = span.c_strand
        int i = 0
        int length = len(segments)
        bint strandprob = False
        bint chromprob = False

    if length > 0:
        seg0 = segments[0]
        if len(chain._segments) == 0:
            my_chrom  = seg0.chrom
            my_strand = seg0.c_strand

        while i < length:
            seg = segments[i]
            if seg.chrom != my_chrom:
                chromprob = True
                break
                
            if seg.c_strand != my_strand:
                strandprob = True
                break

            i += 1

        if strandprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(chain)))
        if chromprob == True:
            raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(chain)))
    return (my_chrom,strand_to_str(my_strand))


cdef class s8(object):

    def __cinit__(self, *segments, **attr):
        cdef:
            str my_chrom, my_strand
            GenomicSegment seg, seg0
            list positions = []
            int length = len(segments)
            int num_segs
            list new_segments
            long total_length
            long c,x
            dict ihash = {}
            array.array my_hash
            array.array pmask

        if "type" not in attr:
            attr["type"] = "exon"
        self.attr = attr
        self._mask_segments = []
        self._segments      = []
        self._inverse_hash    = {}
        self.length = 0
        self.masked_length = 0
        self.spanning_segment = NullSegment
        if length == 0:
            self._position_hash   = array.clone(hash_template,0,False)
            self._position_mask   = array.clone(mask_template,0,False)
        else:
            my_chrom, my_strand = s8_check_segments(self,segments)

            # add new positions
            for seg in segments:
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            new_segments = positionlist_to_segments(my_chrom,my_strand,sorted(set(positions)))
            self._segments = new_segments
            num_segs = len(new_segments)
            
            if num_segs == 0:
                self.spanning_segment = NullSegment 
            elif num_segs == 1:
                self.spanning_segment = new_segments[0]
            elif num_segs >1:
                segs = new_segments
                seg0 = segs[0]
                self.spanning_segment = GenomicSegment(seg0.chrom,
                                                       seg0.start,
                                                       segs[-1].end,
                                                       seg0.strand)

            self._update()
#            length = sum([len(X) for X in segments])
#            my_hash = array.clone(hash_template,length,False)  
#
#            c = 0
#            for seg in new_segments:
#                x = seg.start
#                while x < seg.end:
#                    my_hash[c] = x
#                    ihash[x] = c
#                    c += 1
#                    x += 1
#            
#            self._position_hash = my_hash
#            self._inverse_hash  = ihash
#            self.length = length
#            self.masked_length = length
#            pmask = array.clone(mask_template,length,True)
#            self._position_mask = pmask

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

    cdef void _update(self):
        cdef:
            list segments = self._segments
            int num_segs = len(segments)
            GenomicSegment segment, seg0
            long length = sum([len(X) for X in segments])
            long x, c
            array.array my_hash = array.clone(hash_template,length,False)  
            long [:] my_view = my_hash
            dict ihash = {}

        self.length = length
        #self._position_hash = my_hash
        self._position_mask = array.clone(mask_template,length,True)
        self._inverse_hash  = ihash
        self.c_reset_masks()

#        if num_segs == 0:
#            self.spanning_segment = NullSegment
#        else:
        c = 0
        for segment in segments:
            x = segment.start
            while x < segment.end:
                my_view[c] = x
                ihash[x] = c
                c += 1
                x += 1
        
        if num_segs == 1:
            self.spanning_segment = self._segments[0]
        else:
            segs = self._segments
            seg0 = segs[0]
            self.spanning_segment = GenomicSegment(seg0.chrom,
                                                   seg0.start,
                                                   segs[-1].end,
                                                   seg0.strand)

    cpdef void sort(self):
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
    
    def get_position_list2(self):
        return self._position_hash.aslongs()

    def get_position_list(self):
        return self._position_hash.tolist()

    cdef list c_get_position_list(self, bint copy):
        cdef list positions = self._position_hash.tolist()
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    cdef set c_get_position_set(self):
        return set(self.c_get_position_list(False))
    
    def get_position_set(self):
        return self.c_get_position_set()

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
       
    cdef tuple check_segments(self, tuple segments):
        cdef:
            GenomicSegment seg, seg0
            GenomicSegment span = self.spanning_segment
            str my_chrom  = span.chrom
            Strand my_strand = span.c_strand
            int i = 0
            int length = len(segments)
            bint strandprob = False
            bint chromprob = False

        if length > 0:
            seg0 = segments[0]
            if len(self._segments) == 0:
                my_chrom  = seg0.chrom
                my_strand = seg0.c_strand

            while i < length:
                seg = segments[i]
                if seg.chrom != my_chrom:
                    chromprob = True
                    break
                    
                if seg.c_strand != my_strand:
                    strandprob = True
                    break

                i += 1

            if strandprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches strand (%s) of SegmentChain '%s'" % (seg,my_strand,str(self)))
            if chromprob == True:
                raise ValueError("Incoming GenomicSegment '%s' mismatches chromosome (%s) of SegmentChain '%s'" % (seg,my_chrom,str(self)))
        return (my_chrom,strand_to_str(my_strand))

    cdef void c_add_segments(self, tuple segments):
        cdef:
            str my_chrom, my_strand
            list positions = self.c_get_position_list(False)
            GenomicSegment seg
            int length = len(segments)

        if length > 0:
            my_chrom, my_strand = s8_check_segments(self,segments)
            # add new positions
            for seg in segments: # changed
                positions.extend(range(seg.start,seg.end))
            
            # reset variables
            self._segments = positionlist_to_segments(my_chrom,my_strand,sorted(set(positions)))
            self._update()

    def add_segments(self,*segments):
        if len(segments) > 0:
            if len(self._mask_segments) > 0:
                warnings.warn("Segmentchain: adding segments to %s will reset its masks!",UserWarning)

            self.c_add_segments(segments)
    
    cdef void c_reset_masks(self):
        cdef int [:] pmask = self._position_mask
        pmask [:] = 0
        self._mask_segments = []
        self.masked_length = self.length

    def reset_masks(self):
        self.c_reset_masks()


