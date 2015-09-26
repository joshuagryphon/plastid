# cython: embedsignature=True

import re
import copy
from plastid.genomics.c_common cimport Strand, forward_strand, reverse_strand, \
                                    unstranded, undef_strand

segpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)\(([+-.])\)")
igvpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)")


# compile time constants - __richcmp__ test
DEF LT  = 0
DEF LEQ = 1
DEF EQ  = 2
DEF NEQ = 3
DEF GT  = 4
DEF GEQ = 5


#==============================================================================
# Exported helper functions e.g. for sorting or building |GenomicSegments|
#==============================================================================

cpdef positions_to_segments(str chrom, str strand, object positions):
    """Construct |GenomicSegments| from a chromosome name, a strand, and a list of chromosomal positions

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
        List of |GenomicSegments| covering `positions`
    """
    if isinstance(positions,set):
        positions = sorted(positions)
    else:
        positions = sorted(set(positions))
    return positionlist_to_segments(chrom,strand,positions)

cpdef positionlist_to_segments(str chrom, str strand, list positions):
    """Construct |GenomicSegments| from a chromosome name, a strand, and a list of chromosomal positions

    Parameters
    ----------
    chrom : str
        Chromosome name
       
    strand : str
        Chromosome strand (`'+'`, `'-'`, or `'.'`)
       
    positions : list of unique integers
        Sorted, **end-inclusive** list of positions to include
        in final |GenomicSegment|
           
    Returns
    -------
    list
        List of |GenomicSegments| covering `positions`


     .. warning::
        
        This function is meant to quickly without excessive type conversions.
        So, the elements `positions` must be **UNIQUE** and **SORTED**. If
        they are not, use :func:`positions_to_segments` instead.
    """
    cdef:
        list segments = []
        long start, last_pos, i

    if len(positions) > 0:
        start = positions[0]
        last_pos = positions[0]
        for i in positions:
            if i == start:
                continue
            if i - last_pos == 1:
                last_pos = i
            else:
                segments.append(GenomicSegment(chrom,start,last_pos+1,strand))
                start = i
                last_pos = i
        segments.append(GenomicSegment(chrom,start,last_pos+1,strand))

    return segments

cpdef sort_segments_lexically(GenomicSegment seg):
    """Key function to sort |GenomicSegments| lexically by genomic position,
    by (in order of precedence): chromosome, start, end, strand

    Parameters
    ----------
    seg : |GenomicSegment|

    Returns
    -------
    str
        Chromosome name
        
    int
        Leftmost coordinate of |GenomicSegment|
        
    int
        Rightmost coordinate of |GenomicSegment|
        
    str
        Chromosome strand (`'+'`, `'-'`, or `'.'`)
    """
    return (seg.chrom,seg.start,seg.end,seg.strand)


#==============================================================================
# Helpers
#==============================================================================

cdef void nonecheck(object obj,str place, str valname):
    """Propagate errors if incoming objects are None"""
    if obj is None:
        raise ValueError("%s: value of %s cannot be None" % (place,valname))

cdef str strand_to_str(Strand strand):
    """Convert enum Strand to str representation"""
    if strand == forward_strand:
        return "+"
    elif strand == reverse_strand:
        return "-"
    elif strand == unstranded:
        return "."
    elif strand == undef_strand:
        return "strand undefined"
    else:
        raise ValueError("strand_to_str: Strand must be forward (%s), reverse (%s), or unstranded(%s). Got '%s'" % (forward_strand,
            reverse_strand, unstranded, strand))

cdef Strand str_to_strand(str val):
    """Convert str representation of strand to enum"""
    if val == "+":
        return forward_strand
    elif val == "-":
        return reverse_strand
    elif val == ".":
        return unstranded
    elif val == "\x00":
        return undef_strand
    else:
        raise ValueError("Strand must be '+', '-', '.', or '\\x00' (undefined)")


#==============================================================================
# Classes
#==============================================================================

cdef class GenomicSegment:
    """A continuous segment of the genome, defined by a chromosome name,
    a start coordinate, and end coordinate, and a strand. Building block
    for a |SegmentChain| or a |Transcript|.

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
   

    Examples
    --------

    |GenomicSegments| sort lexically by chromosome, start position, end position,
    and finally strand::
        
        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrB",0,10,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrA",75,100,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrA",55,75,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrA",50,150,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrA",50,100,"-")
        True


    They also provide a few convenience methods for containment or overlap. To be
    contained, a segment must be on the same chromosome and strand as its container,
    and its coordinates must be within or equal to its endpoints

        >>> GenomicSegment("chrA",50,100,"+") in GenomicSegment("chrA",25,100,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") in GenomicSegment("chrA",50,100,"+")
        True
       
        >>> GenomicSegment("chrA",50,100,"+") in GenomicSegment("chrA",25,100,"-")
        False

        >>> GenomicSegment("chrA",50,100,"+") in GenomicSegment("chrA",75,200,"+")
        False


    Similarly, to overlap, |GenomicSegments| must be on the same strand and 
    chromosome.
     
    See also
    --------
    SegmentChain
        Base class for genomic features, built from multiple |GenomicSegments|

    Transcript
        Subclass of |SegmentChain| that adds convenience methods for manipulating
        coding regions, UTRs, et c
    """
    def __cinit__(self, str chrom, long start, long end, str strand):
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
        if end < start:
            raise ValueError("GenomicSegment: start coordinate (%s) must be >= end (%s)." % (start,end))
        my_chrom = chrom
        self.chrom  = my_chrom
        self.start  = start
        self.end    = end
        self.c_strand = str_to_strand(strand)
    
    def __reduce__(self): # enables pickling
        return (GenomicSegment,(self.chrom,self.start,self.end,self.strand))

    def __repr__(self):
        return "<%s %s:%s-%s strand='%s'>" % ("GenomicSegment",
                                              self.chrom,
                                              self.start,
                                              self.end,
                                              self.strand)
    
    def __str__(self):
        return "%s:%s-%s(%s)" % (self.chrom, self.start, self.end, self.strand)
    
    @staticmethod
    def from_str(str inp):
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
        cdef:
            long start, end
            str chrom, strand, s_start, s_end
        chrom, s_start, s_end, strand = segpat.search(inp).groups()
        start = long(s_start)
        end   = long(s_end) 
        return GenomicSegment(chrom,start,end,strand)
    
    def __copy__(self):
        cdef str newchrom, newstrand
        cdef long newstart, newend
        newchrom = copy.copy(self.chrom)
        newstrand = self.strand
        newstart = self.start
        newend = self.end
        return GenomicSegment(newchrom,newstart,newend,newstrand)

    def __deepcopy__(self,memo):
        return self.__copy__()

    def __len__(self):
        """Return length, in nucleotides, of |GenomicSegment|"""
        cdef long len
        len = self.end - self.start
        return len

    def __richcmp__(self,GenomicSegment other, int cmptype):
        if other is None or not isinstance(other,GenomicSegment):
            return False

        return self._cmp_helper(other,cmptype)

    # TODO: suspend type check via decorator, since we have already done this
    cpdef bint _cmp_helper(self,GenomicSegment other,int cmptype):
        nonecheck(other,"GenomicSegment eq/neq","other")
        schrom = self.chrom
        ochrom = other.chrom
        if cmptype == EQ:
            return schrom         == ochrom and\
                   self.c_strand  == other.c_strand and\
                   self.start     == other.start and\
                   self.end       == other.end
        elif cmptype == NEQ:
            return self._cmp_helper(other,EQ) == False
        elif cmptype == LT:
            if schrom < ochrom:
                return True
            elif schrom > ochrom:
                return False
            elif schrom == ochrom:
                sstart = self.start
                ostart = other.start
                if sstart < ostart:
                    return True
                elif sstart > ostart:
                    return False
                elif sstart == ostart:
                    send = self.end
                    oend = other.end
                    if send < oend:
                        return True
                    elif send > oend:
                        return False
                    elif send == oend:
                        return self.c_strand < other.c_strand
        elif cmptype == GT:
            return other._cmp_helper(self,LT) # (other < self)
        elif cmptype == LEQ:
            return self._cmp_helper(other,LT) or self._cmp_helper(other,EQ) # (self == other) or (self < other)
        elif cmptype == GEQ:
            return other._cmp_helper(self,LT) or self._cmp_helper(other,EQ) # (self == other) or (other < self)
        else:
            raise AttributeError("Comparison operation not defined")
    
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
        nonecheck(other,"GenomicSegment.__eq__","other")
        return self.contains(other)

    cpdef bint contains(self,GenomicSegment other):
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
        nonecheck(other,"GenomicSegment.contains","other")
        return self.chrom == other.chrom and\
               self.c_strand == other.c_strand and\
               (other.start >= self.start and other.end <= self.end and other.end >= other.start)
         
    cpdef bint overlaps(self,GenomicSegment other):
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
        nonecheck(other,"GenomicSegment.overlaps","other")
        if self.chrom == other.chrom and self.c_strand == other.c_strand:
            if (self.start >= other.start and self.start < other.end) or\
               (other.start >= self.start and other.start < self.end):
                   return True

        return False

    cpdef str as_igv_str(self):
        """Format as an IGV location string"""
        return "%s:%s-%s" % (self.chrom, self.start+1, self.end+1)
    
    @staticmethod
    def from_igv_str(str loc_str, str strand="."):
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
        cdef:
            str chrom, s_start, s_end
            long start, end

        chrom,s_start,s_end = igvpat.search(loc_str).groups()
        start = long(s_start) - 1
        end   = long(s_end) - 1
        return GenomicSegment(chrom,start,end,strand)
    
    property start:
        """Zero-indexed (Pythonic) start coordinate of |GenomicSegment|"""
        def __get__(self):
            return self.start

    property end:
        """Zero-indexed, half-open (Pythonic) end coordinate of |GenomicSegment|"""
        def __get__(self):
            return self.end

    property chrom:
        """Chromosome where |GenomicSegment| resides"""
        def __get__(self):
            return self.chrom

    property strand:
        """Strand of |GenomicSegment|:

          - '+' for forward / Watson strand
          - '-' for reverse / Crick strand
          - '.' for unstranded / both strands
        """
        def __get__(self):
            return strand_to_str(self.c_strand)

    property c_strand:
        def __get__(self):
            return self.c_strand

# Exported object
NullSegment = GenomicSegment("NullChromosome",0,0,"\x00")
