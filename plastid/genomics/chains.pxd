from plastid.genomics.c_roitools cimport GenomicSegment, Strand
from plastid.genomics.c_common cimport ExBool, bool_exception
from cpython cimport array
import numpy
import array
cimport numpy


cdef class s1(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        long [:] _position_hash
        int  [:] _position_mask
        dict attr, _inverse_hash
    
    cpdef void _update(self)
    cdef void _get_position_hash(self)
    cpdef void sort(self)
    cdef numpy.ndarray c_get_position_list(self, bint)
    cpdef set get_position_set(self)
    cdef tuple check_segments(self, tuple)
    cdef c_add_segments(self, tuple)
    cpdef void reset_masks(self)

cdef class s2(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        long [:] _position_hash
        int  [:] _position_mask
        dict attr, _inverse_hash

    cpdef void _update(self)
    cdef void _get_position_hash(self)
    cpdef void sort(self)
    cdef numpy.ndarray c_get_position_list(self, bint)
    cpdef set get_position_set(self)
    cdef tuple check_segments(self, tuple)
    cdef c_add_segments(self, tuple)
    cpdef void reset_masks(self)

cdef class s3(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        long [:] _position_hash
        int  [:] _position_mask
        dict attr, _inverse_hash

    cpdef void _update(self)
    cdef void _get_position_hash(self)
    cpdef void sort(self)
    cdef numpy.ndarray c_get_position_list(self, bint)
    cpdef set get_position_set(self)
    cdef tuple check_segments(self, tuple)
    cdef c_add_segments(self, tuple)
    cpdef void reset_masks(self)

cdef tuple check_segments(s2,tuple)
cdef void add_segments(s2,tuple)

cdef tuple s3_check_segments(s3,tuple)
cdef void s3_setup(s3,tuple)

# fastest init: 10.4-10.9us
cdef class s4(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        long [:] _position_hash
        int  [:] _position_mask
        dict attr, _inverse_hash

    cpdef void _update(self)
    cdef void _get_position_hash(self)
    cpdef void sort(self)
    cdef numpy.ndarray c_get_position_list(self, bint)
    cpdef set get_position_set(self)
    cdef tuple check_segments(self, tuple)
    cdef void c_add_segments(self, tuple)
    cpdef void reset_masks(self)

cdef class s5(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        long [:] _position_hash
        int  [:] _position_mask
        dict attr, _inverse_hash

    cpdef void _update(self)
    cdef void _get_position_hash(self)
    cpdef void sort(self)
    cdef numpy.ndarray c_get_position_list(self, bint)
    cpdef set get_position_set(self)
    cdef tuple check_segments(self, tuple)
    cdef c_add_segments(self, tuple)
    cpdef void reset_masks(self)


cdef class s6(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        long [:] _position_hash
        int  [:] _position_mask
        dict attr, _inverse_hash

    cdef void _update(self)
    cpdef void sort(self)
    cdef numpy.ndarray c_get_position_list(self, bint)
    cdef set c_get_position_set(self)
    cdef tuple check_segments(self, tuple)
    cdef void c_add_segments(self, tuple)
    cdef void c_reset_masks(self)


cdef tuple s6_check_segments(s6, tuple)

# matches s6 for init empty
# 10% faster for add_segments. much faster if using arrays instead of memoryviews
# slower for init with segments. can be faster with arrays instead of memoryviews?
cdef class s7(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        #long [:] _position_hash # faster as array.array
        #int [:] _position_mask  # faster as array.array
        array.array _position_hash, _position_mask
        dict attr, _inverse_hash

    cdef void _update(self)
    cpdef void sort(self)
    cdef list c_get_position_list(self, bint)
    cdef set c_get_position_set(self)
    cdef tuple check_segments(self, tuple)
    cdef void c_add_segments(self, tuple)
    cdef void c_reset_masks(self)

cdef class s8(object):
    cdef:
        list _segments, _mask_segments
        GenomicSegment spanning_segment
        long length, masked_length
        #long [:] _position_hash # faster as array.array
        #int [:] _position_mask  # faster as array.array
        array.array _position_hash, _position_mask
        dict attr, _inverse_hash

    cdef void _update(self)
    cpdef void sort(self)
    cdef list c_get_position_list(self, bint)
    cdef set c_get_position_set(self)
    cdef tuple check_segments(self, tuple)
    cdef void c_add_segments(self, tuple)
    cdef void c_reset_masks(self)


