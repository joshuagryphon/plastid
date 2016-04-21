cdef class CenterMapFactory:
    cdef unsigned int nibble

cdef class FivePrimeMapFactory:
    cdef int offset

cdef class ThreePrimeMapFactory:
    cdef int offset

cdef class VariableFivePrimeMapFactory:
    cdef int [10000] forward_offsets
    cdef int [10000] reverse_offsets

cdef class StratifiedVariableFivePrimeMapFactory(VariableFivePrimeMapFactory):
    cdef int min_length, max_length, _numlengths
     

cdef class SizeFilterFactory:
    cdef int min_, max_