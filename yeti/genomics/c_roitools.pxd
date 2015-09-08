cdef enum Strand:
   forward_strand, reverse_strand, unstranded, problem

cdef class GenomicSegment(object):
    cdef char* chrom
    cdef long start
    cdef long end
    cdef Strand strand
    cpdef bint contains(self,GenomicSegment)
    cpdef bint overlaps(self,GenomicSegment)
    cpdef str as_igv_str(self)
    cpdef str get_name(self)


cdef positionlist_to_segments(str,str,list)
