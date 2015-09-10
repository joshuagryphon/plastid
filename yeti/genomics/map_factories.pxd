cdef class CenterMapFactory(object):
    cdef unsigned int nibble

cdef class FivePrimeMapFactory(object):
    cdef int offset

cdef class ThreePrimeMapFactory(object):
    cdef int offset

cdef class VariableFivePrimeMapFactory(object):
    cdef int [10000] forward_offsets
    cdef int [10000] reverse_offsets

cdef class SizeFilterFactory(object):
    cdef float min_, max_

