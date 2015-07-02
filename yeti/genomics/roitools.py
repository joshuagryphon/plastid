#!/usr/bin/env python
"""This module defines object types that describe features or regions of interest in a genome or contig.


Important classes
-----------------
|GenomicSegment|
    A fundamental unit of a feature, similar to :py:class:`HTSeq.GenomicInterval`.
    |GenomicSegment| describes a single region of a genome, and is specified by
    by a chromosome name, a start coordinate, an end coordinate, and a strand.
    |GenomicSegment| contain provide no feature annotation data, and are used
    primarily to construct |SegmentChain| or |Transcript| objects, which do
    provide feature annotation data. |GenomicSegment| implements various methods
    to test equality to, overlap with, and containment of other |GenomicSegment|
    objects.

|SegmentChain|
    Base class for genomic features with annotation data. |SegmentChains|
    can contain zero or more |GenomicSegments|, and therefore can model
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
"""
__date__ = "2011-09-01"
__author__ = "joshua"
import numpy
import re
import copy

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from yeti.util.services.decorators import skipdoc
from yeti.util.services.colors import get_str_from_rgb255, get_rgb255_from_str
from yeti.readers.gff_tokens import make_GFF3_tokens, \
                                    make_GTF2_tokens



igvpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)")
segpat  = re.compile(r"([^:]*):([0-9]+)-([0-9]+)\(([+-.])\)")
ivcpat = re.compile(r"([^:]*):([^(]+)\(([+-.])\)")


#===============================================================================
# io functions
#===============================================================================

@skipdoc
def _get_attr(dtmp):
    ltmp=[]
    for k,v in sorted(dtmp.items()):
        ltmp.append("%s='%s'" % (k,v))
    
    return ",".join(ltmp)

@skipdoc
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

@skipdoc
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

def positionlist_to_segments(chrom,strand,positions):
    """Construct |GenomicSegment| from a chromosome name, a strand, and a list of chromosomal positions

    Parameters
    ----------
    chrom : str
        Chromosome name
       
    strand : str
        Chromosome strand (`'+'`, `'-'`, or `'.'`)
       
    positions : list of integers
        **End-inclusive** list, tuple, or set of positions to include
        in final |GenomicSegment|
           
    Returns
    -------
    list
        List of |GenomicSegment| objects covering `positions`
    """
    positions = sorted(list(set(positions)))
    intervals = []
    if len(positions) > 0:
        last_pos = positions[0]
        start = positions[0]
        for i in sorted(positions[1:]):
            if i - last_pos == 1:
                last_pos = i
            else:
                intervals.append(GenomicSegment(chrom,start,last_pos+1,strand))
                start = i
                last_pos = i
        intervals.append(GenomicSegment(chrom,start,last_pos+1,strand))
    return intervals


#===============================================================================
# sorting functions
#===============================================================================

def sort_segments_lexically(seg):
    """Key function to sort a list of |GenomicSegment| lexically
    by genomic position, by (in order of precedence): chromosome, start, end, strand

    Parameters
    ----------
    seg : |GenomicSegment|

    Returns
    -------
    str : Chromosome name
    int : Leftmost coordinate of |GenomicSegment|
    int : Rightmost coordinate of |GenomicSegment|
    str : Chromosome strand (`'+'`, `'-'`, or `'.'`)
    """
    return (seg.chrom,seg.start,seg.end,seg.strand)

def sort_segmentchains_lexically(segchain):
    """Key function to sort a list of features lexically by genomic
    position, by (in order of precedence): chromosome, start, end, strand,
    then by length (in nucleotides) and name should any of the other
    criteria fail

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
    """Key function to sort |SegmentChain| by name"""
    return segchain.get_name()


#===============================================================================
# classes that handle genome annotation features
#===============================================================================

class GenomicSegment(object):
    """A continuous segment of the genome, defined by a chromosome name,
    a start coordinate, and end coordinate, and a strand.

    Attributes
    ----------
    chrom : str
        Name of chromosome

    start : int
        0-indexed, left most position of segment

    end : int
        0-indexed, half-open right most position of segment

    strand : str
        Chromosome strand (`'+'`, `'-'`, or `'.'`)
    
    
    See also
    --------
    SegmentChain
        Base class for genomic features, built from multiple |GenomicSegments|
    """
    
    __slots__ = ("chrom","start","end","strand")
    
    def __init__(self,chrom,start,end,strand):
        """Create a |GenomicSegment|
        
        Parameters
        ----------
        chrom : str
            Chromosome name
        
        start : int
            0-indexed, leftmost coordinate of feature
        
        end : int
            0-indexed, half-open rightmost coordinate of feature
            Must be >= `start`
        
        strand : str
            Chromosome strand (`'+'`, `'-'`, or `'.'`)
        """
        assert strand in ("+","-",".")
        assert end >= start
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
    
    def __repr__(self):
        return "<%s %s:%s-%s strand='%s'>" % (self.__class__.__name__,
                                                  self.chrom,
                                                  self.start,
                                                  self.end,
                                                  self.strand)
    
    def __str__(self):
        return "%s:%s-%s(%s)" % (self.chrom, self.start, self.end, self.strand)
    
    @staticmethod
    def from_str(inp):
        """Construct a |GenomicSegment| from its ``str()`` representation
        
        Parameters
        ----------
        inp : str
            String representation of |GenomicSegment| as `chrom:start-end(strand)`
            where `start` and `end` are in 0-indexed, half-open coordinates
        
        Returns
        -------
        |GenomicSegment|
        """
        chrom, start, end, strand = segpat.search(inp).groups()
        start = int(start)
        end   = int(end) 
        return GenomicSegment(chrom,start,end,strand)
    
    def __len__(self):
        """Return length, in nucleotides, of |GenomicSegment|"""
        return self.end - self.start
    
    def __eq__(self,other):
        """Test whether `self` and `other` cover the exact same set of positions,
        on the same chromosome and strand.
        
        Parameters
        ----------
        other : |GenomicSegment|
            Query segment
        
        Returns
        -------
        bool
        """
        return isinstance(other,GenomicSegment) and\
               self.chrom  == other.chrom and\
               self.strand == other.strand and\
               self.start  == other.start and\
               self.end    == other.end
    
    def __ne__(self,other):
        return False if self == other else True
    
    def __contains__(self,other):
        """Test whether this segment contains `other`, where containment is
        defined as all positions in `other` being present in self, when both
        `self` and `other` share the same chromosome and strand.
           
        Parameters
        ----------
        other : |GenomicSegment|
            Query segment
        
        Returns
        -------
        bool
        """
        my_range = set(range(self.start,self.end+1))
        return self.chrom == other.chrom and\
               self.strand == other.strand and\
               other.start in my_range and\
               other.end in my_range
    
    def contains(self,other):
        """Test whether this segment contains `other`, where containment is
        defined as all positions in `other` being present in self, when both
        `self` and `other` share the same chromosome and strand.
           
        Parameters
        ----------
        other : |GenomicSegment|
            Query segment
        
        Returns
        -------
        bool
        """
        return other in self
        
    def overlaps(self,other):
        """Test whether this segment overlaps `other`, where overlap is defined
        as sharing: a chromosome, a strand, and a subset of coordinates.
           
        Parameters
        ----------
        other : |GenomicSegment|
            Query segment
        
        Returns
        -------
        bool
        """
        if self.chrom == other.chrom and self.strand == other.strand:
            my_range  = set(range(self.start,self.end))
            its_range = set(range(other.start,other.end))
            if len(my_range & its_range) > 0:
                    return True
            else:
                return False
        else:
            return False
    
    def as_igv_str(self):
        """Format as an IGV location string"""
        return "%s:%s-%s" % (self.chrom, self.start+1, self.end+1)
    
    @staticmethod
    def from_igv_str(loc_str,strand="."):
        """Construct |GenomicSegment| from IGV location string
        
        Parameters
        ----------
        igvloc : str
            IGV location string, in format `'chromosome:start-end'`,
            where `start` and `end` are 1-indexed and half-open
            
        strand : str
            The chromosome strand (`'+'`, `'-'`, or `'.'`)
            
        Returns
        -------
        |GenomicSegment|
        """
        chrom,start,end = igvpat.search(loc_str).groups()
        start = int(start) - 1
        end   = int(end) - 1
        return GenomicSegment(chrom,start,end,strand)
    
    def get_name(self):
        """Alias for :meth:`GenomicSegment.__str__`"""
        return str(self)


#===============================================================================
# higher-order classes that handle multi-feature structures, like transcripts
# or alignments
#===============================================================================

class SegmentChain(object):
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
        self.spanning_segment = None  #interval spanning entire SegmentChain
        self._segments      = []
        self._mask_segments = []      # list<GenomicSegment> of masked positions
        self._position_hash = self._get_position_hash()
        self.strand = None
        self.chrom  = None
        self.attr = attr
        
        if "type" not in attr:
            self.attr["type"] = "exon"
        
        self.add_segments(*segments)

    def _update(self):
        """Sort component |GenomicSegments| within the |SegmentChain|,
        and maintain synchrony of position hashes
        """
        self.sort()
        self._position_hash = self._get_position_hash()
        if len(self) > 0:
            self.spanning_segment = GenomicSegment(self.chrom,
                                      self[0].start,
                                      self[-1].end,
                                      self.strand)

    def sort(self):
        """Sort component segments by ascending 5' chromosomal coordinate"""
        self._segments.sort(key=sort_segments_lexically)
        self._mask_segments.sort(key=sort_segments_lexically)

    def __repr__(self):
        sout = "<%s intervals=%s" % (self.__class__.__name__, len(self))
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
    
    def __setitem__(self,key,val):
        self._segments[key] = val
            
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
    
    def __contains__(self,other):
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
        if isinstance(other,GenomicSegment):
            tmp = SegmentChain(other)
            return tmp in self
        elif isinstance(other,SegmentChain):        
            if len(self) == 0 or len(other) == 0:
                return False
            elif other.get_length() > self.get_length():
                return False
            elif self.chrom != other.chrom:
                return False
            elif self.strand != other.strand:
                return False
            elif len(other) == 1:
                for segment in self:
                    if segment.contains(other[0]):
                        return True
                return False            
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
                    return False
                
                if found is True:
                    return ojuncs == myjuncs[mystart:mystart+len(ojuncs)]
                else:
                    return False
        else:
            raise TypeError("Other is not a GenomicSegment or SegmentChain")
        return False
    
    def __eq__(self,other):
        """Test whether `self` and `other` are equal. Equality is defined as
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
            other = SegmentChain(GenomicSegment)
            
        if len(self) == 0 or len(other) == 0:
            return False
        else:
            return self.chrom == other.chrom and\
                   self.strand == other.strand and\
                   self.get_position_set() == other.get_position_set()
        
#         if len(self) != len(other):
#             return False
#         elif len(self) == 0:
#             return False
#         elif len(self) > 1:
#             return self in other and other in self
#         else:
#             return self[0].start == other[0].start and\
#                    self[0].end == other[0].end and\
#                    self.strand == other.strand and\
#                    self.chrom == other.chrom
                   
    def __ne__(self,other):
        """Defines inequality for two SegmentChains as complement of :meth:`SegmentChain.__eq__`

        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        Returns
        -------
        bool
        """
        return False if self == other else True

    def shares_segments_with(self,other):
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
            if self.chrom != other.chrom or self.strand != other.strand:
                return []
            else:
                shared = []
                for segment in other:
                    if segment in self._segments:
                        shared.append(segment)
                return shared
        else:
            raise TypeError("Other is not a GenomicSegment or SegmentChain")
    
    def unstranded_overlaps(self,other):
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
        if isinstance(other,GenomicSegment):
            tmp = SegmentChain(other)
            return self.unstranded_overlaps(tmp)
        elif isinstance(other,SegmentChain):
            if self.chrom != other.chrom:
                return False
            else:
                my_pos = self.get_position_set()
                o_pos  = other.get_position_set()
            return len(my_pos & o_pos) > 0 
        else:
            raise TypeError("Other is not a GenomicSegment or SegmentChain")
        
    def overlaps(self,other):
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
        if isinstance(other,GenomicSegment):
            tmp = SegmentChain(other)
            return self.overlaps(tmp)
        return self.strand == other.strand and self.unstranded_overlaps(other) 
    
    def antisense_overlaps(self,other):
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
        if isinstance(other,GenomicSegment):
            tmp = SegmentChain(other)
            return self.antisense_overlaps(tmp)
        return self.strand != other.strand and self.unstranded_overlaps(other) 

    def covers(self,other):
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
        if isinstance(other,GenomicSegment):
            return self.covers(SegmentChain(other))
        elif len(self) == 0 or len(other) == 0:
            return False
        elif isinstance(other,SegmentChain):
            return self.strand == other.strand and\
                   self.chrom  == other.chrom and\
                   other.get_position_set() & self.get_position_set() == other.get_position_set()
        else:
            raise TypeError("Other is not a GenomicSegment or SegmentChain")
    
    def get_antisense(self):
        """Returns an |SegmentChain| antisense to `self`, with empty `attr` dict.
        
        Returns
        -------
        SegmentChain
            |SegmentChain| antisense to `self`
        """
        new_strand = "+" if self.strand == "-" else "-" if self.strand == "+" else "."
        new_segments = [GenomicSegment(X.chrom,X.start,X.end,new_strand) for X in self]
        return SegmentChain(*tuple(new_segments))

    def get_position_list(self):
        """Retrieve a sorted list of genomic coordinates included in this |SegmentChain|
        
        Returns
        -------
        list<int>
        """
        return sorted(list(self.get_position_set()))
    
    def get_position_set(self):
        """Retrieve a set of genomic coordinates included in this |SegmentChain|
        
        Returns
        -------
        set<int>
        """
        positions = set(self._position_hash.keys())
        return positions
        
    def _get_position_hash(self):
        """Creates a hash mapping relevant genomic positions to |SegmentChain| coordinates
        
        Returns
        -------
        dict
            Dictionary mapping each position in the |SegmentChain| to its 
            chromosomal coordinate
        """
        self.sort()
        my_hash = {}
        c = 0
        for segment in self:
            for x in range(segment.start,segment.end):
                my_hash[x] = c
                c += 1
        return my_hash
    
    def get_valid_position_set(self):
        """Returns a set of genomic coordinates which are not covered by any masks

        Returns
        -------
        set<int>
        """
        position_set = self.get_position_set()
        masked = []
        for segment in self._mask_segments:
            masked.extend(range(segment.start,segment.end))
        
        masked = set(masked)
        return position_set - masked
    
    def get_name(self):
        """Returns the name of this |SegmentChain|, first searching through
        ``self.attr`` for the keys ``ID``, ``Name``, and ``name``.
        If one is not found, a name is generated one from the locations in the
        collection.
        
        Returns
        -------
        str
            Returns in order of preference, ID from ``self.attr``, 'Name' from
            ``self.attr``, "name" from ``self.attr`` or str(self) 
        """
        name = self.attr.get("ID",
               self.attr.get("Name",
               self.attr.get("name",
                             str(self))))
        return name
    
    def get_gene(self):
        """Return name of gene associated with |SegmentChain|, if any, 
        by searching through ``self.attr`` for the keys ``gene_id`` and ``Parent``.
        If one is not found, a generated gene name for the SegmentChain is 
        made from :py:meth:`get_name`.

        Returns
        -------
        str
            Returns in order of preference, 'gene_id' from ``self.attr``, 
            'Parent' from ``self.attr`` or ``"gene_%s" % self.get_name()``
        """
        gene = self.attr.get("gene_id",
               self.attr.get("Parent",
               "gene_%s" % self.get_name()))
        if isinstance(gene,list):
            gene = ",".join(sorted(gene))
            
        return gene
    
    def get_length(self):
        """Return total length, in nucleotides, of this |SegmentChain|
        
        Returns
        -------
        int
        """
        return sum([len(X) for X in self])

    def get_valid_length(self):
        """Return the total length, in nucleotides, of all unmasked
        positions in this |SegmentChain|
        
        Returns
        -------
        int
        """
        return len(self.get_valid_position_set())

    def add_segment(self,segment):
        """Add a |GenomicSegment| to the |SegmentChain|. If there are
        already intervals in the collection, the incoming segments must be 
        on the same strand and chromosome as all others present.
        
        Parameters
        ----------
        segment : GenomicSegment
            Interval to add to |SegmentChain|
        """
        self.add_segments(segment)
        
    def add_segments(self,*segments):
        """Add 1 or more |GenomicSegment| to the |SegmentChain|. If there are
        already intervals in the collection, the incoming segments must be 
        on the same strand and chromosome as all others present.

        Parameters
        ----------
        segments : |GenomicSegment|
            One or more |GenomicSegment| to add to |SegmentChain|
        """
        if len(segments) > 0:
            strands = set([X.strand for X in segments])
            chroms  = set([X.chrom for X in segments])
            assert len(strands) == 1
            assert len(chroms)  == 1
            if len(self) > 0:
                assert segments[0].chrom  == self.chrom
                assert segments[0].strand == self.strand
            else:
                self.chrom  = segments[0].chrom
                self.strand = segments[0].strand
            
            # add new positions
            positions = self.get_position_set()
            for segment in segments:
                positions |= set(range(segment.start,segment.end))
             
            self._segments = positionlist_to_segments(self.chrom,self.strand,positions)
        self._update()
        
    def add_masks(self,*mask_segments):
        """Adds one or more |GenomicSegment| to the collection of masks. If there are
        already segments in the collection, the incoming segments must be 
        on the same strand and chromosome as all others present.

        Parameters
        ----------
        mask_segments : |GenomicSegment|
            One or more segments, in genomic coordinates, covering positions to exclude
			from return values of :py:meth:`get_valid_position_set()`,
			:py:meth:`get_valid_counts()`, or :py:meth:`get_valid_length()`
		
		See also
		--------
		SegmentChain.get_masks
		
        SegmentChain.reset_masks
        """
        if len(mask_segments) > 0:
            strands = set([X.strand for X in mask_segments])
            chroms  = set([X.chrom for X in mask_segments])
            assert len(strands) == 1
            assert len(chroms)  == 1
            if len(self) > 0:
                assert mask_segments[0].chrom  == self.chrom
                assert mask_segments[0].strand == self.strand
            else:
                self.chrom  = mask_segments[0].chrom
                self.strand = mask_segments[0].strand
            
            # add new positions to any existing masks
            positions = set()
            for segment in list(mask_segments) + self._mask_segments:
                positions |= set(range(segment.start,segment.end))
             
            # trim away non-overlapping masks
            positions &= self.get_position_set()
            
            # regenerate list of ivs from positions, in case some were doubly-listed
            self._mask_segments = positionlist_to_segments(self.chrom,self.strand,positions)
        self._update()
    
    def get_masks(self):
        """Return masked positions as a list of |GenomicSegment|
        
        Returns
        -------
        list
            list of |GenomicSegment| s representing masked positions
        
        See also
        --------
        SegmentChain.add_masks
        
        SegmentChain.reset_masks
        """
        return self._mask_segments
    
    def get_masks_as_IVC(self):
        """Return masked positions as an |SegmentChain|
        
        Returns
        -------
        |SegmentChain|
            Masked positions
        
        See also
        --------
        SegmentChain.add_masks"thickend",

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
            
    def get_junctions(self):
        """Returns a list of |GenomicSegment| representing spaces
        between the |GenomicSegment| in self._segments. In the case of a transcript,
        these would represent introns. In the case of an alignment, these
        would represent gaps in the query compared to the reference.
        
        Returns
        -------
        list<|GenomicSegment|>
            List of |GenomicSegment| covering spaces between the intervals in ``self``
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
    
    def as_gff3(self,feature_type=None,escape=True,excludes=[]):
        """Format a length-1 |SegmentChain| as a line of `GFF3`_ output.
        
        Because `GFF3`_ files permit many schemas of parent-child hierarchy,
        and in order to reduce confusion and overhead, attempts to export
        a multi-interval |SegmentChain| will raise an :py:obj:`AttributeError`.
        
        Instead, users may export the individual features from which the
        multi-interval |SegmentChain| was constructed, or construct features
        for them, setting *ID*, *Parent*, and *type* attributes following
        their own conventions.
        
        
        Columns of `GFF3`_ are as follows:
        
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
        
        More info can be found at:
            - `GFF3 file format specification <http://www.sequenceontology.org/gff3.shtml>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_

         
        Parameters
        ----------
        feature_type : str
            If not None, overrides the *type* attribute of ``self.attr``
        
        escape : bool, optional
            Escape tokens in column 9 of `GFF3`_ output (Default: *True*)
        
        excludes : list, optional
            List of attribute key names to exclude from column 9
            (Default: *[]*)
        
        Returns
        -------
        str
            Line of `GFF3`_-formatted text
        
            
        Raises
        -----
        AttributeError
            if the SegmentChain has multiple intervals
        """ 
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
                           "type",]
        
        for segment in self:
            ltmp2 = self._get_8_gff_columns(segment,feature_type) +\
                   [make_GFF3_tokens(gff_attr,
                                     excludes=always_excluded+excludes,
                                     escape=True)]

        return "\t".join(ltmp2) + "\n"
        
    
    def as_gtf(self,feature_type=None,escape=True,excludes=[]):
        """Format |SegmentChain| as a block of GTF2 output.
        
        The ``frame`` or ``phase`` attribute (GTF2 column 8) is valid only for "CDS"
        features, and, if not present in ``self.attr``, is calculated assuming
        the |SegmentChain| contains the entire coding region. If the |SegmentChain|
        contains multiple intervals, the ``frame`` or ``phase`` attribute will
        *always* be recalculated.
        
        All attributes in ``self.attr``, except those created upon import,
        will be propagated to all of the features that are generated.
        
        Columns of `GTF2`_ are as follows:
        
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
        
        More info can be found at:
            - `GTF2 file format specification <http://mblab.wustl.edu/GTF22.html>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_

         
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
        The `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_ requires
        that attributes ``gene_id`` and ``transcript_id`` be defined. If these
        are not present in ``self.attr``, their values will be guessed 
        following the rules in :py:meth:`SegmentChain.get_gene` and 
        :py:meth:`SegmentChain.get_name`, respectively.
        
        To save memory, only the attributes shared by all of the individual
        sub-features (e.g. exons) that were used to assemble this |SegmentChain|
        have been stored in ``self.attr``. This means that upon re-export to GTF,
        these sub-features will be lacking any attributes that were specific
        to them individually. Formally, this is compliant with the 
        `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_, which states
        explicitly that only the attributes ``gene_id`` and ``transcript_id``
        are supported.
        """
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
                           "color"]
        
        for segment in self:
            ltmp2 = self._get_8_gff_columns(segment,feature_type) +\
                   [make_GTF2_tokens(gtf_attr,
                                     excludes=always_excluded+excludes,
                                     escape=escape)]
            
            ltmp1.append("\t".join(ltmp2))            

        return "\n".join(ltmp1) + "\n"
    
    def _get_8_gff_columns(self,segment,feature_type):
        """Format columns 1-8 of GFF/GTF2/GFF3 files. 
        
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
        phase = "."
        if feature_type == "CDS":
            # use phase/frame if known for length-1 features
            # called "phase" in GFF3 conventions; "frame" in GTF2
            if len(self) == 1 and ("phase" in self.attr or "frame" in self.attr):
                phase = self.attr.get("phase",self.attr.get("frame"))
            # otherwise calculate
            else:
                segment_start = self.get_segmentchain_coordinate(segment.chrom,segment.start,segment.strand,stranded=False)
                phase = (3 - (segment_start % 3)) % 3
        
        ltmp = [segment.chrom,
                self.attr.get("source","."),
                feature_type,
                segment.start + 1,
                segment.end + 1 - 1,
                str(self.attr.get("score",".")),
                segment.strand,
                phase]
        
        return [str(X) for X in ltmp]
        
        
    
    def as_bed(self,thickstart=None,thickend=None,as_int=True,color=None):
        """Format |SegmentChain| as a string of BED12 output.

        BED12 columns are as follows (see
        the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        for more details).

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


        Parameters
        ----------
        as_int : bool, optional
            Force "score" to integer (Default: True)
    
        thickstart : int or None, optional
            If not None, overrides the genome coordinate to start
            plotting as thick in genome browser found in ``self.attr["thickstart"]``
    
        thickend : int or None, optional
            If not None, overrides the genome coordinate to stop
            plotting as thick in genome browser found in ``self.attr["thickend"]``
    
        color : str or None, optional
            Color represented as RGB hex string.
            If not none, overrides the color in ``self.attr["color"]``
    
    
        Returns
        -------
        str
            Line of BED12-formatted text
        """
        if len(self) > 0:
            score = self.attr.get("score",0)
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
            
            thickstart = self.attr.get("thickstart",self[0].start) if thickstart is None else thickstart
            thickend   = self.attr.get("thickend",self[0].start)   if thickend   is None else thickend
            
            ltmp = [self[0].chrom,
                    self[0].start,
                    self[-1].end,
                    self.get_name(),
                    score,
                    self[0].strand,
                    thickstart,
                    thickend,
                    color,
                    len(self),
                    ",".join([str(len(X)) for X in self]) + ",",
                    ",".join([str(X.start - self[0].start) for X in self]) + ","
                   ]            
            return "\t".join([str(X) for X in ltmp]) + "\n"
        else:
            # SegmentChain with no intervals
            return ""
    
    def as_psl(self):
        """Formats |SegmentChain| as PSL (blat) output.
        
        Notes
        -----
        This will raise an :py:class:`AttributeError` unless the following
        keys are present and defined in ``self.attr``, corresponding to the
        columns of a PSL file:
        
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
        
        
        Returns
        -------
        str
            PSL-representation of BLAT alignment

        
        Raises
        ------
        AttributeError
            If not all of the attributes listed above are defined
        """
        ltmp = []
        try:
            ltmp.append(self.attr["match_length"])
            ltmp.append(self.attr["mismatches"])
            ltmp.append(self.attr["rep_matches"])
            ltmp.append(self.attr["N"])
            ltmp.append(self.attr["query_gap_count"])
            ltmp.append(self.attr["query_gap_bases"])
            ltmp.append(self.attr["target_gap_count"])
            ltmp.append(self.attr["target_gap_bases"])
            ltmp.append(self.attr["strand"])
            ltmp.append(self.attr["query_name"])
            ltmp.append(self.attr["query_length"])
            ltmp.append(self.attr["query_start"])
            ltmp.append(self.attr["query_end"])
            ltmp.append(self.attr["target_name"])
            ltmp.append(self.attr["target_length"])
            ltmp.append(self.attr["target_start"])
            ltmp.append(self.attr["target_end"])
            ltmp.append(len(self))
    
            block_sizes = ",".join([str(len(X)) for X in self]) + ","
            q_starts = ",".join([str(X) for X in self.attr["q_starts"]]) + ","
            t_starts = ",".join([str(X) for X in self.attr["t_starts"]]) + ","
            
            ltmp.append(block_sizes)
            ltmp.append(q_starts)
            ltmp.append(t_starts)
            return "\t".join(str(X) for X in ltmp) + "\n"
        except KeyError:
            raise AttributeError("SegmentChains only support PSL output if all PSL attributes are defined in self.attr: match_length, mismatches, rep_matches, N, query_gap_count, query_gap_bases, strand, query_length, query_start, query_end, target_name, target_length, target_start, target_end")
        
    def get_segmentchain_coordinate(self,chrom,genomic_x,strand,stranded=True):
        """Finds the |SegmentChain| coordinate corresponding to a genomic position
        
        Parameters
        ----------
        chrom : str
            Chromosome name
            
        genomic_x : int
            x coordinate, in genomic space
            
        strand : str
            Chromosome strand ('+', '-', or '.')
            
        stranded : bool, optional
            if ``True``, coordinates are given in stranded space
            (what you'd expect for a transcript)
        
        
        Returns
        -------
        int
            Position in |SegmentChain|
            
        Raises
        ------
        KeyError : if position outside bounds of |SegmentChain|
        """
        assert chrom  == self.chrom
        assert strand == self.strand
        if self._position_hash == {}:
            self._position_hash = self._get_position_hash()
        if stranded is True and self.strand == "-":
            return self.get_length() - self._position_hash[genomic_x] - 1
        else:
            return self._position_hash[genomic_x]
    
    def get_genomic_coordinate(self,x,stranded=True):
        """Finds genomic coordinate corresponding to position ``x`` in this |SegmentChain|
        
        Parameters
        ----------
        x : int
            position of interest, relative to |SegmentChain|
            
        stranded : bool
            If true and the |SegmentChain| is on the minus strand,
            ``x`` will be given relative to the end coordinate of the 
            |SegmentChain| rather than the start (default: True)
        
                             
        Returns
        -------
        str 
            Chromosome name
        
        int
            Genomic cordinate corresponding to ``x``
        
        str
            Chromosome strand ('+', '-', or '.')
        
        
        Raises
        ------
        IndexError
            if ``x`` is outside the bounds of the |SegmentChain|
        """
        positions = self.get_position_list()
        if stranded is True and self.strand == "-":
            positions = positions[::-1]
        return self.chrom, positions[x], self.strand
    
    def get_subchain(self,start,end,stranded=True):
        """Retrieves a sub-|SegmentChain| corresponding a range of positions
        specified in coordinates relative this |SegmentChain|. ``self.attr`` is
        copied to the child SegmentChain.
        
        Parameters
        ----------
        start : int
            position of interest in SegmentChain coordinates, 0-indexed
            
        end : int
            position of interest in SegmentChain coordinates, 0-indexed 
            and half-open
            
        stranded : bool
            If ``True`` and the SegmentChain is on the minus strand,
            `x` will be given relative to the end coordinate of the 
            SegmentChain rather than the start (default: ``True``)
            
                          
        Returns
        -------
        |SegmentChain|
            covering positions `start` to `end` of parent |SegmentChain|
        
        
        Raises
        ------
        IndexError
            if `start` or `end` is outside the bounds of the |SegmentChain|
        """
        if start is None or end is None:
            raise IndexError()
        positions = self.get_position_list()
        if stranded is True and self.strand == "-":
            positions = positions[::-1]
        positions = positions[start:end]
        ivs = positionlist_to_segments(self.chrom,self.strand,set(positions))
        return SegmentChain(*tuple(ivs),**copy.deepcopy(self.attr))

    def get_counts(self,gnd,stranded=True):
        """Returns list of counts or values at each position in SegmentChain
        
        Parameters
        ----------
        gnd : Subclass of AbstractGenomeArray
            GenomeArray from which to fetch counts
            
        stranded : bool
            If true and the SegmentChain is on the minus strand,
            count order will be reversed relative to genome
            (default: ``True``)
            
            
        Returns
        -------
        list<Number> 
        """
        ltmp = []
        for iv in self:
            ltmp.extend(gnd[iv])
        if self.strand == "-" and stranded is True:
            ltmp = ltmp[::-1]
        return ltmp
    
    def get_valid_counts(self,gnd,stranded=True):
        """Returns counts of |SegmentChain| as a masked array, in transcript
		coordinates. Positions masked by :py:meth:`SegmentChain.add_mask`
		will be masked
        
        Parameters
        ----------
        gnd : non-abstract subclass of |AbstractGenomeArray|
            GenomeArray from which to fetch counts
            
        stranded : bool
            If true and the |SegmentChain| is on the minus strand,
            count order will be reversed relative to genome
            (default: ``True``)
            
            
        Returns
        -------
		:py:class:`numpy.ma.masked_array`
        """
        ltmp = []
        for iv in self:
            ltmp.extend(gnd[iv])
        if self.strand == "-" and stranded is True:
            ltmp = ltmp[::-1]
        
        atmp = numpy.ma.masked_invalid(ltmp)
        atmp.mask = True
        
        valid_positions = [self.get_segmentchain_coordinate(self.spanning_segment.chrom,X,
                                                            self.spanning_segment.strand) 
                           for X in self.get_valid_position_set()]
        for x in valid_positions:
            atmp.mask[x] = False
        return atmp
    
    def get_sequence(self,genome,stranded=True):
        """Returns spliced genomic sequence of |SegmentChain| as a string
        
        Parameters
        ----------
        genome : dict
            Dictionary mapping chromosome names to sequences.
            Sequences may be strings, string-like, or :py:class:`Bio.Seq.SeqRecord` objects
       
        stranded : bool
            If true and the |SegmentChain| is on the minus strand,
            sequence will be reverse-complemented (default: True)
            
            
        Returns
        -------
        str
            Nucleotide sequence of the |SegmentChain| evaluated over ``genome``
        """
        chromseq = genome[self.spanning_segment.chrom]
        if not isinstance(chromseq,SeqRecord):
            chromseq = SeqRecord(Seq(chromseq,generic_dna))
            
        ltmp = [chromseq[X.start:X.end] for X in self]
        stmp = "".join([str(X.seq) for X in ltmp])
        if self.strand == "-"  and stranded is True:
            stmp = str(Seq(stmp,generic_dna).reverse_complement())
            
        return stmp
    
    def get_fasta(self,genome,stranded=True):
        """Formats sequence of SegmentChain as FASTA output
        
        Parameters
        ----------
        genome : dict
            Dictionary mapping chromosome names to sequences.
            Sequences may be strings, string-like, or :py:class:`Bio.Seq.SeqRecord` objects
       
        stranded : bool
            If true and the |SegmentChain| is on the minus strand,
            sequence will be reverse-complemented (default: True)

            
        Returns
        -------
        str
            FASTA-formatted seuqence of |SegmentChain| evaluated over ``genome``
        """
        return ">%s\n%s\n" % (self.get_name(),self.get_sequence(genome,stranded=stranded))

    @staticmethod
    def from_str(inp):
        """Factory method to create an SegmentChains from a lines from a string
		formatted as by :py:meth:`SegmentChain.str`:
           
            chrom:start-end^start-end(strand)
           
        where '^' indicates a splice junction between regions specified
        by `start` and `end` and `strand` is '+', '-', or '.'. Coordinates are
        0-indexed and half-open


        Parameters
        ----------
        inp : str
			String formatted in manner of :py:meth:`SegmentChain.str`
          
          
        Returns
        -------
        |SegmentChain|
        """
        if inp in ("na","nan","None:(None)","None","none",None) or isinstance(inp,float) and numpy.isnan(inp):
            return SegmentChain()
        else:
            chrom,middle,strand = ivcpat.search(inp).groups()
            ivs = []
            for piece in middle.split("^"):
                start,end = piece.split("-")
                start = int(start)
                end = int(end)
                ivs.append(GenomicSegment(chrom,start,end,strand))
            return SegmentChain(*tuple(ivs))
        
    @staticmethod
    def from_bed(line):
        """Create |SegmentChain| from a line from a BED file.
        
        See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        for more details.

        Parameters
        ----------
        line
            Line from a BED file, containing 4 or more columns


        Returns
        -------
        |SegmentChain|
        """
        frags = []
        items = line.strip("\n").split("\t")
        chrom         = items[0]
        chrom_start   = int(items[1])
        chrom_end     = int(items[2])
        strand = "." if len(items) < 6 else items[5]
    
        default_id  = "%s:%s-%s(%s)" % (chrom,chrom_start,chrom_end,strand)
    
        # dict mapping optional bed column to tuple of (Name,default value)
        # these values are used if any optional columns 4-12 are ommited
        bed_columns = { 3 :  ("ID",         default_id,    str),
                        4 :  ("score",      numpy.nan,     float),
                        #5 :  ("strand",     strand),
                        6 :  ("thickstart", None,          int),
                        7 :  ("thickend",   None,          int),
                        8 :  ("color",      "0,0,0",       str),
                        9 :  ("blocks",     "1",             int),
                        10 : ("blocksizes", str(chrom_end - chrom_start),str),
                        11 : ("blockstarts","0",             str),
                      }
    
        # set attr defaults in case we're dealing with BED4-BED9 format
        attr = { KEY : VAL for KEY,VAL,_ in bed_columns.values() }
    
        # populate attr with real values from columns that are present
        for i, tup in sorted(bed_columns.items()):
            if len(items) > i:
                key     = tup[0]
                default = tup[1]
                func    = tup[2]
                try:
                    attr[key] = func(items[i])
                except ValueError:
                    attr[key] = default
            else:
                break
    
        # convert color to hex string
        try:
            attr["color"] = get_str_from_rgb255(tuple([int(X) for X in attr["color"].split(",")]))
        except ValueError:
            attr["color"] = "#000000"
    
        # sanity check on thickstart and thickend
        if attr["thickstart"] == attr["thickend"]: # if coding region is 0 length, RNA is non-coding
            attr["thickstart"] = attr["thickend"] = None
        elif any([attr["thickstart"] is None, attr["thickend"] is None]):
            attr["thickstart"] = attr["thickend"] = None
        elif attr["thickstart"] < 0 or attr["thickend"] < 0:
            attr["thickstart"] = attr["thickend"] = None
        
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
    
        return SegmentChain(*tuple(frags),**attr)
    
    @staticmethod
    def from_psl(psl_line):
        """Create an |SegmentChain| from a line from a PSL (BLAT) file
        
        Parameters
        ----------
        psl_line : str
            Line from a PSL file
        
        See the `PSL spec <http://pombe.nci.nih.gov/genome/goldenPath/help/blatSpec.html>`_
        
        Returns
        -------
        |SegmentChain|
        """
        items = psl_line.strip().split("\t")        
        attr = {}
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
        
        ivs = []
        for t_start, block_size in zip(t_starts,block_sizes):
            iv = GenomicSegment(attr["target_name"],
                                 t_start,
                                 t_start + block_size,
                                 attr["strand"])
            ivs.append(iv)
        
        return SegmentChain(*ivs,**attr)        
    
class Transcript(SegmentChain):
    """Subclass of |SegmentChain| specifically modeling transcripts.
    In addition to coordinate-conversion, count fetching, sequence fetching,
    and various other methods inherited from |SegmentChain|, |Transcript|
    provides convenience methods for fetching sub-regions containing only
    CDS features, 5' UTRs, and 3' UTRs.

    Intervals are sorted from lowest to greatest starting coordinate on their
    reference sequence, regardless of strand. Iteration over the SegmentChain
    will yield intervals from left-to-right in the genome.
    
    Tests for equality, inequality, containment, and overlapping  are provided.

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
        """Create a Transcript
        
        Parameters
        ----------
        attr : dict
            dictionary of attributes

        ivs : |GenomicSegment|
            Any number >=0 of |GenomicSegment|
        
        attr : dict
            Any relevant metadata

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
            If provided, a transcript_id used for GTF2 export.
            Otherwise, generated from genomic coordinates.
        
        attr["gene_id"] : str
            If provided, a gene_id used for GTF2 export
            Otherwise, generated from genomic coordinates.
        """
        self._segments   = []
        if "type" not in attr:
            attr["type"] = "mRNA"
            
        self.cds_genome_start = attr.get("cds_genome_start",None)
        self.cds_genome_end   = attr.get("cds_genome_end",None)
        SegmentChain.__init__(self,*ivs,**attr)
        self._update()
 
    def _update(self):
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
        """Returns the name of this |SegmentChain|, first searching through
        ``self.attr`` for the keys ``transcript_id``, ``ID``, ``Name``, and ``name``.
        If one is not found, a name is generated one from the locations in the
        collection.
        
        Returns
        -------
        str
            Returns in order of preference, ID from ``self.attr``, 'Name' from
            ``self.attr``, "name" from ``self.attr`` or str(self) 
        """
        name = self.attr.get("transcript_id",
               self.attr.get("ID",
               self.attr.get("Name",
               self.attr.get("name",
                             str(self)))))
        return name
    
    def get_cds(self):
        """Retrieve sub-SegmentChain covering CDS, including the stop codon.
        If no coding region, returns an empty |SegmentChain|.
        
        The following attributes are passed from ``self.attr`` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`
        
        Returns
        -------
        |SegmentChain|
            CDS region of |Transcript| if present, otherwise empty |SegmentChain|
        """
        my_segmentchain = SegmentChain()
        if self.cds_genome_start is not None and self.cds_genome_end is not None:
            my_segmentchain = Transcript(*self.get_subchain(self.cds_start,self.cds_end,stranded=True),
                                cds_genome_start=self.cds_genome_start,
                                cds_genome_end=self.cds_genome_end,
                                transcript_id=self.get_name(),
                                type="CDS",gene_id=self.get_gene())
            
        return my_segmentchain   
    
    def get_utr5(self):
        """Retrieve sub-|SegmentChain| covering 5'UTR of |Transcript|
        If no coding region, returns an empty |SegmentChain|

        The following attributes are passed from ``self.attr`` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`
        
        Returns
        -------
        |SegmentChain|
            5' UTR region of |Transcript| if present, otherwise empty |SegmentChain|
        """
        my_segmentchain = SegmentChain()
        if self.cds_genome_start is not None and self.cds_genome_end is not None:
            my_segmentchain = self.get_subchain(0,self.cds_start,stranded=True)
            
            my_segmentchain.attr["type"] = "5UTR"
            my_segmentchain.attr["gene_id"] = self.get_gene()
            my_segmentchain.attr["transcript_id"] = self.get_gene()
        return my_segmentchain   
    
    def get_utr3(self):
        """Retrieve sub-SegmentChain covering 3'UTR of |Transcript|, excluding
        the stop codon. If no coding region, returns an empty |SegmentChain|
        
        The following attributes are passed from ``self.attr`` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`

        Returns
        -------
        |SegmentChain|
            3' UTR region of |Transcript| if present, otherwise empty |SegmentChain|
        """
        my_segmentchain = SegmentChain()
        if self.cds_genome_start is not None and self.cds_genome_end is not None:
            my_segmentchain = self.get_subchain(self.cds_end,self.get_length(),
                                                   stranded=True)
            my_segmentchain.attr["type"] = "3UTR"
            my_segmentchain.attr["gene_id"] = self.get_gene()
            my_segmentchain.attr["transcript_id"] = self.get_gene()
            
        return my_segmentchain   

    def as_gtf(self,feature_type="exon",escape=True):
        """Format Transcript as a GTF2 block. |GenomicSegment| are formatted
        as GTF2 "exon" features. Coding regions, if peresent, are formatted
        as GTF2 "CDS" features. Stop codons are excluded in the "CDS" features,
        per the GTF2 specification, and exported separately.

        All attributes from ``self.attr`` are propagated to the exon and CDS
        features that are generated.

        Columns of GTF2 are as follows:
        
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
        
        More info can be found at:
            - `GTF2 file format specification <http://mblab.wustl.edu/GTF22.html>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
         
         
        Parameters
        ----------
        feature_type : str
            If not None, overrides the "type" attribute of self.attr
        
        escape : bool, optional
            URL escape tokens in column 9 of GTF output (Default: True)
        
        
        Returns
        -------
        str
            Block of GTF2-formatted text


        Notes
        -----
        The `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_ requires
        that attributes ``gene_id`` and ``transcript_id`` be defined. If these
        are not present in ``self.attr``, their values will be guessed 
        following the rules in :py:meth:`SegmentChain.get_gene` and 
        :py:meth:`SegmentChain.get_name`, respectively.
        
        To save memory, only the attributes shared by all of the individual
        sub-features (e.g. exons) that were used to assemble this |SegmentChain|
        have been stored in ``self.attr``. This means that upon re-export to GTF,
        these sub-features will be lacking any attributes that were specific
        to them individually. Formally, this is compliant with the 
        `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_, which states
        explicitly that only the attributes ``gene_id`` and ``transcript_id``
        are supported.
        """
        stmp  = SegmentChain.as_gtf(self,feature_type=feature_type,escape=escape)
        cds_ivc_temp = self.get_cds()
        if len(cds_ivc_temp) > 0:
            child_ivc_attr  = copy.deepcopy(self.attr)
            child_ivc_attr.pop("type")
            cds_positions = cds_ivc_temp.get_position_list()
            
            # remove stop codons from CDS, per GTF2 spec
            if self.spanning_segment.strand == "+":
                cds_positions = cds_positions[:-3]
            else:
                cds_positions = cds_positions[3:]
            cds_ivc = SegmentChain(*positionlist_to_segments(self.spanning_segment.chrom,
                                                        self.spanning_segment.strand,
                                                        cds_positions),
                                   type="CDS",**child_ivc_attr)
            stmp += cds_ivc.as_gtf(feature_type="CDS",escape=escape)
            
            start_codon_ivc = cds_ivc_temp.get_subchain(0, 3)
            start_codon_ivc.attr.update(child_ivc_attr)
            stmp += start_codon_ivc.as_gtf(feature_type="start_codon",escape=escape)
    
            stop_codon_ivc = cds_ivc_temp.get_subchain(cds_ivc_temp.get_length()-3, cds_ivc_temp.get_length())
            stop_codon_ivc.attr.update(child_ivc_attr)
            stmp += stop_codon_ivc.as_gtf(feature_type="stop_codon",escape=escape)

        return stmp

    def as_gff3(self,escape=True,excludes=[],rna_type="mRNA"):
        """Format a |Transcript| as a block of `GFF3`_ output, following
        the schema set out in the `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`
        
        The |Transcript| will be formatted according to the following rules:
        
          1. A feature of type ``rna_type`` will be created, with *Parent* attribute
             set to the value of :meth:`Transcript.get_gene`, and *ID* attribute
             set to :meth:`Transcript.get_name`
        
          2. For each |GenomicSegment| in the |Transcript|, a child feature of type
             *exon* will be created. The *Parent* attribute of these features
             will be set to the value of :meth:`Transcript.get_name`. These will
             have unique IDs generated from :meth:`Transcript.get_name`.

          3. If the |Transcript| is coding (i.e. has none-*None* value for
             ``self.cds_genome_start`` and ``self.cds_genome_end``), child features
             of type *five_prime_UTR*, *CDS*, and *three_prime_UTR* will be created,
             with *Parent* attributes set to :meth:`Transcript.get_name`. These will
             have unique IDs generated from :meth:`Transcript.get_name`.
        
        
        Parameters
        ----------
        escape : bool, optional
            Escape tokens in column 9 of `GFF3`_ output (Default: *True*)
        
        excludes : list, optional
            List of attribute key names to exclude from column 9
            (Default: *[]*)
        
        rna_type : str, optional
            Feature type to export RNA as (e.g. *"tRNA"*, *"noncoding_RNA"*,
            et c. Default: *"mRNA"*)

        
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
            `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`,
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
        gene_id       = self.get_gene()
        transcript_id = self.get_name()
        ltmp = []
        
        child_attr = copy.deepcopy(self.attr)
        keys_to_pop = ("ID",)
        for k in keys_to_pop:
            if k in child_attr:
                child_attr.pop(k)

        # mRNA feature
        feature = SegmentChain(self.spanning_segment,ID=transcript_id,Parent=gene_id,type=rna_type)
        ltmp.append(feature.as_gff3(excludes=excludes,escape=escape))

        # child features
        child_attr["Parent"] = transcript_id

        # exon feature
        child_attr["type"] = "exon"
        for n,iv in enumerate(self):
            my_id   = "%s:exon:%s" % (transcript_id,n)
            child_attr["ID"]   = my_id
            feature = SegmentChain(iv,**child_attr)
            ltmp.append(feature.as_gff3(excludes=excludes,escape=escape))
        
        # CDS & UTRs
        if self.cds_genome_start is not None:
            parts = [("five_prime_UTR", self.get_utr5()),
                     ("CDS",            self.get_cds()),
                     ("three_prime_UTR",self.get_utr3()),
                    ]
            for ftype, ivc in parts:
                child_attr["type"]   = ftype
                child_attr["Parent"] = transcript_id
                for n,iv in enumerate(ivc):
                    my_id   = "%s:%s:%s" % (transcript_id,ftype,n)
                    child_attr["ID"]   = my_id
                    feature = SegmentChain(iv,**child_attr)
                    ltmp.append(feature.as_gff3(excludes=excludes,escape=escape))

        return "".join(ltmp)
    
    def as_bed(self,as_int=True,color=None):
        """Formats |Transcript| as a BED12 line, assigning CDS boundaries 
        to the thickstart and thickend columns from ``self.attr``

        BED12 columns are as follows (see
        the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        for more details).

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


        Parameters
        ----------
        as_int : bool, optional
            Force "score" to integer (Default: True)
    
        color : str or None, optional
            Color represented as RGB hex string.
            If not none, overrides the color in `self.attr["color"]`
    
    
        Returns
        -------
        str
            Line of BED12-formatted text
        """
        return SegmentChain.as_bed(self,
                                   thickstart=self.cds_genome_start,
                                   thickend=self.cds_genome_end,
                                   as_int=as_int,
                                   color=color)

    @staticmethod
    def from_bed(line):
        """Factory method to create a |Transcript| from a BED line,
        assuming that the thickstart and thickend columns specify CDS boundaries
    
    	See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
    	for more details.

        Parameters
        ----------
        line
            Line from a BED file with at least 4 columns
    
    
        Returns
        -------
        |Transcript|
        """
        ivc = SegmentChain.from_bed(line)
        segments = ivc._segments
        attr = ivc.attr
        attr.pop("type") # default type for SegmentChain is "exon". We want to use "mRNA"
        attr["cds_genome_start"] = attr["thickstart"]
        attr["cds_genome_end"]   = attr["thickend"]
        attr.pop("thickstart")
        attr.pop("thickend")
        transcript = Transcript(*segments,**attr)
    
        return transcript
    
    @staticmethod
    def from_psl(psl_line):
        ivc = SegmentChain.from_psl(psl_line)
        segments = ivc._segments
        attr = ivc.attr
        attr["cds_genome_start"] = None
        attr["cds_genome_end"] = None
        return Transcript(*segments,**attr)
        
