cdef enum Strand:
    undef_strand   = 0
    forward_strand = 1   # can use bitwise operations to detect overlap
    reverse_strand = 2
    unstranded     = 3

cdef enum ExBool:
    bool_exception = -1
    false     = 0
    true      = 1

