from plastid.genomics.c_common cimport Strand, undef_strand, ExBool, bool_exception
from cpython cimport array
import numpy
import array
cimport numpy


cdef str strand_to_str(Strand)
cdef Strand str_to_strand(str) #except undef_strand

cpdef list positions_to_segments(str,str,object)
cpdef list positionlist_to_segments(str,str,list)
cpdef Transcript add_three_for_stop_codon(Transcript)
cpdef list merge_segments(list)

cdef class GenomicSegment(object):                      # 20 bytes + str
    cdef str chrom
    cdef long start                                     # 8 bytes
    cdef long end                                       # 8 bytes
    cdef Strand c_strand                                # 4 bytes if int
    cpdef bint contains(self,GenomicSegment)
    cpdef bint overlaps(self,GenomicSegment)
    cpdef bint _cmp_helper(self,GenomicSegment,int)
    cpdef str as_igv_str(self)

cdef class SegmentChain(object):                       # >= 836 bytes
    cdef:
        list _segments, _mask_segments                 # 72 bytes/each -> 144 + contents
        readonly GenomicSegment spanning_segment       # 20 bytes + str
        readonly long length                           # 8 bytes
        readonly long masked_length                    # 8 bytes

        array.array _position_hash                     # 56 + contents
        array.array _position_mask # really should be bint, but we can't use it there
        array.array _endmap
        public dict attr                               # 280 bytes + contents
        #dict _inverse_hash                             # 280 bytes + contents

    # maintenance of chain internals
    cdef bint c_add_segments(self,tuple) except False
    cdef bint _set_segments(self,list) except False
    cdef array.array _get_position_hash(self)
    cdef bint _set_masks(self,list) except False
    cdef void c_reset_masks(self)
    #cdef dict _get_inverse_hash(self)

    # comparison operators
    cdef  list c_shares_segments_with(self, SegmentChain)
    cdef  ExBool c_covers(self ,SegmentChain) except bool_exception
    cdef  ExBool c_antisense_overlaps(self, SegmentChain) except bool_exception
    cdef  ExBool c_unstranded_overlaps(self, SegmentChain) except bool_exception 
    cdef  ExBool c_overlaps(self, SegmentChain) except bool_exception
    cdef  ExBool c_contains(self, SegmentChain) except bool_exception

    # position info - c impl
    cdef list c_get_position_list(self)
    cdef set c_get_position_set(self)
    cdef numpy.ndarray c_get_position_array(self,bint)
    cpdef set get_masked_position_set(self)

    # subchains/coordinates
    cdef long c_get_genomic_coordinate(self, long, bint) except -1
    cdef long c_get_segmentchain_coordinate(self, long, bint) except -1
    cdef SegmentChain c_get_subchain(self,long, long, bint)
    cpdef SegmentChain get_antisense(self)

 
cdef class Transcript(SegmentChain):
    cdef:
        object cds_genome_start, cds_genome_end, cds_start, cds_end

    cdef bint _update_cds(self) except False
    cdef bint _update_from_cds_start(self) except False
    cdef bint _update_from_cds_end(self) except False
    cdef bint _update_from_cds_genome_start(self) except False
    cdef bint _update_from_cds_genome_end(self) except False

cdef GenomicSegment NullSegment
