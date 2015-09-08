cdef enum Strand:
   forward_strand, reverse_strand, unstranded, problem

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


cpdef sort_segments_lexically(GenomicSegment)
cpdef positions_to_segments(str,str,object)
cpdef positionlist_to_segments(str,str,list)
