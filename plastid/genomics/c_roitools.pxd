from plastid.genomics.c_common cimport Strand, undef_strand

cdef class GenomicSegment(object):
    cdef str chrom
    cdef long start
    cdef long end
    cdef Strand c_strand
    cpdef bint contains(self,GenomicSegment)
    cpdef bint overlaps(self,GenomicSegment)
    cpdef bint _cmp_helper(self,GenomicSegment,int)
    cpdef str as_igv_str(self)

cdef str strand_to_str(Strand)
cdef Strand str_to_strand(str) #except undef_strand

cpdef sort_segments_lexically(GenomicSegment)
cpdef positions_to_segments(str,str,object)
cpdef positionlist_to_segments(str,str,list)

