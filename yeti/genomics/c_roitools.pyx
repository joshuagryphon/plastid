import re
import copy

segpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)\(([+-.])\)")
igvpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)")


cdef positionlist_to_segments(str chrom, str strandstr, list positions):
    positions = sorted(list(set(positions)))
    cdef list segments = []
    cdef long start = positions[0]
    cdef long last_pos = positions[0]
    cdef long i
    cdef Strand strand = str_to_strand(strandstr)
    if len(positions) > 0:
        start = positions[0]
        last_pos = positions[0]
        for i in positions:
            if i == start:
                continue
            if i - last_pos == 1:
                last_pos = 1
            else:
                segments.append(GenomicSegment(chrom,start,last_pos+1,strand))
                start = i
                last_pos = i
        segments.append(GenomicSegment(chrom,start,last_pos+1,strand))

    return segments

cdef void nonecheck(object obj,str place, str valname):
    if obj is None:
        raise ValueError("%s: value of %s cannot be None" % (place,valname))

cdef str strand_to_str(Strand strand):
    """Convert enum Strand to str representation"""
    if strand == forward_strand:
        return "+"
    elif strand == reverse_strand:
        return "-"
    else:
        return "."

cdef Strand str_to_strand(str val) except problem:
    """Convert str representation of strand to enum"""
    if val == "+":
        return forward_strand
    elif val == "-":
        return reverse_strand
    elif val == ".":
        return unstranded
    else:
        raise ValueError("trand must be '+', '-', or '.'")


cdef class GenomicSegment:

    def __cinit__(self, str chrom, long start, long end, str strand):
        if end < start:
            raise ValueError("GenomicSegment: start coordinate must be >= end.")
        my_chrom = chrom
        self.chrom  = my_chrom
        self.start  = start
        self.end    = end
        self.strand = str_to_strand(strand)
    
    def __repr__(self):
        return "<%s %s:%s-%s strand='%s'>" % ("GenomicSegment",
                                              self.chrom,
                                              self.start,
                                              self.end,
                                              strand_to_str(self.strand))
    
    def __str__(self):
        return "%s:%s-%s(%s)" % (self.chrom, self.start, self.end, 
                                 strand_to_str(self.strand))
    
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
        newstrand = strand_to_str(self.strand)
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
        return self._cmp_helper(other,cmptype)

    cpdef bint _cmp_helper(self,GenomicSegment other,int cmptype):
        nonecheck(other,"GenomicSegment eq/neq","other")
        if other is None:
            return False
        if cmptype == 2: # equal 
            return isinstance(other,GenomicSegment) and\
                   self.chrom  == other.chrom and\
                   self.strand == other.strand and\
                   self.start  == other.start and\
                   self.end    == other.end
        elif cmptype == 3: #inequal
            return not self._cmp_helper(other,2)
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
               self.strand == other.strand and\
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
        if self.chrom == other.chrom and self.strand == other.strand:
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
    
    cpdef str get_name(self):
        """Alias for :meth:`GenomicSegment.__str__`"""
        return str(self)

    property start:
        """Zero-indexed (Pythonic) start coordinate of |GenomicSegment|"""
        def __get__(self):
            return self.start
        def __set__(self, val):
            self.start = <long?>val

    property end:
        """Zero-indexed, half-open (Pythonic) end coordinate of |GenomicSegment|"""
        def __get__(self):
            return self.end
        def __set__(self, val):
            self.end = <long?>val

    property chrom:
        """Chromosome where |GenomicSegment| resides"""
        def __get__(self):
            return self.chrom
        def __set__(self, str val):
            self.chrom = val

    property strand:
        """Strand of |GenomicSegment|:

          - '+' for forward / Watson strand
          - '-' for reverse / Crick strand
          - '.' for unstranded / both strands
        """
        def __get__(self):
            return strand_to_str(self.strand)
        def __set__(self, str val):
            self.val = str_to_strand(val)




