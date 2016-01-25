cdef str strand_to_str(Strand strand):
    """Convert enum Strand to str representation"""
    if strand == forward_strand:
        return "+"
    elif strand == reverse_strand:
        return "-"
    elif strand == unstranded:
        return "."
    elif strand == undef_strand:
        return "strand undefined"
    else:
        raise ValueError("strand_to_str: Strand must be forward (%s), reverse (%s), or unstranded(%s). Got '%s'." % (forward_strand,
            reverse_strand, unstranded, strand))

cdef Strand str_to_strand(str val):
    """Convert str representation of strand to enum"""
    if val == "+":
        return forward_strand
    elif val == "-":
        return reverse_strand
    elif val == ".":
        return unstranded
    elif val == "\x00":
        return undef_strand
    else:
        raise ValueError("Strand must be '+', '-', '.', or '\\x00' (undefined). Got '%s'." % val)

cdef class _GeneratorWrapper(object):
    """Hack to allow generators created in Cython to have REPRs instead of erroring upon inspection"""

    def __cinit__(self,generator,message):
        self.generator = generator
        self.message = "<wrapped generator object '%s'>" % message
        
    def __str__(self):
        return self.message
    
    def __repr__(self):
        return self.message
    
    def __getattr__(self,attr):
        return getattr(self.generator,attr)
