from yeti.genomics.c_roitools cimport GenomicSegment, Strand
from yeti.genomics.c_common cimport ExBool, bool_exception
from cpython import array
import array

cdef class SegmentChain(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        long [:] _position_hash
        int  [:] _position_mask # really should be bint, but we can't use it there
        dict attr, _inverse_hash

    cdef void _update(self)
    #cdef array.array _get_position_hash(self)
    cdef long [:]  _get_position_hash(self)
    cdef dict _update_inverse_hash(self)

    cdef tuple check_segments(self, tuple)
    cdef c_add_segments(self, tuple)

    cpdef void sort(self)

    cdef  list c_shares_segments_with(self, SegmentChain)
    cdef  ExBool c_covers(self ,SegmentChain) except bool_exception
    cdef  ExBool c_unstranded_overlaps(self, SegmentChain) except bool_exception 
    cdef  ExBool c_overlaps(self, SegmentChain) except bool_exception
    cdef  ExBool c_antisense_overlaps(self, SegmentChain) except bool_exception
    cdef  ExBool c_contains(self, SegmentChain) except bool_exception

    cpdef  ExBool c_richcmp(self,SegmentChain,int) except bool_exception

    cpdef SegmentChain get_antisense(self)
    cpdef list get_position_list(self)
    cpdef set get_position_set(self)

    cpdef long get_length(self)
    cdef long c_get_genomic_coordinate(self, long, bint) except -1
    cdef long c_get_segmentchain_coordinate(self, long, bint) except -1
    cdef SegmentChain c_get_subchain(self,long, long, bint)


cdef class Transcript(SegmentChain):
    cdef void _update(self)

