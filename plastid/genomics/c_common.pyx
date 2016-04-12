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

cdef Strand str_to_strand(str val) except error_strand:
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
    """Wrapper class to prevent `repr()` of Cython generators from erroring in interactive sessions.
    All attributes and behaviors are delegated to the wrapped generator, except the attribute
    `generator`.
    
    Attributes
    ----------
    generator : generator
        Wrapped generator object
    """
        
    def __cinit__(self,generator,desc):
        """Wrapper class to prevent `repr()` of Cython generators from erroring in interactive sessions
        
        Parameters
        ----------
        generator : Generator object to wrap
        
        desc : str
            Short description of generator contents, used in generating
            the  `repr` text
        """
        self.generator = generator
        self.desc = "<wrapped generator object '%s'>" % desc
        
    def __str__(self):
        return self.desc
    
    def __repr__(self):
        return self.desc
    
    def __getattr__(self,attr):
        if attr == "generator":
            return self.generator
        else:
            return getattr(self.generator,attr)
    
    def __iter__(self):
        return iter(self.generator)
     
    def __next__(self):
        return next(self.generator)
     
    def next(self):
        return next(self.generator)