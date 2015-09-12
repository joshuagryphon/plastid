cdef enum Strand:
    forward_strand = 0
    unstranded     = 2
    reverse_strand = 1
    undef_strand   = -1

cdef class GenomicSegment(object):
    cdef str chrom
    cdef long start
    cdef long end
    cdef Strand strand
    cpdef bint contains(self,GenomicSegment)
    cpdef bint overlaps(self,GenomicSegment)
    cpdef bint _cmp_helper(self,GenomicSegment,int)
    cpdef str as_igv_str(self)
    cpdef str get_name(self)


cdef str strand_to_str(Strand)
cdef Strand str_to_strand(str) except undef_strand

cpdef sort_segments_lexically(GenomicSegment)
cpdef positions_to_segments(str,str,object)
cpdef positionlist_to_segments(str,str,list)

