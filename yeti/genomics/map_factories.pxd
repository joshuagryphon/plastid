from yeti.genomics.c_common import fusedreal

cdef class CenterMapFactory:
    cdef unsigned int nibble

cdef class FivePrimeMapFactory:
    cdef int offset

cdef class ThreePrimeMapFactory:
    cdef int offset

cdef class VariableFivePrimeMapFactory:
    cdef int [10000] forward_offsets
    cdef int [10000] reverse_offsets

cdef class SizeFilterFactory:
    cdef float min_, max_

