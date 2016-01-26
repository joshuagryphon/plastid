cdef enum Strand:
    error_strand   = 666
    undef_strand   = 0
    forward_strand = 1   # can use bitwise operations to detect overlap
    reverse_strand = 2
    unstranded     = 3

cdef enum ExBool:
    bool_exception = -1
    false     = 0
    true      = 1

cdef str strand_to_str(Strand)
cdef Strand str_to_strand(str) except error_strand

cdef class _GeneratorWrapper(object):
    cdef:
        object generator
        str    desc