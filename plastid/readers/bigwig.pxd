"""Fast reader for `BigWig`_ files, implemented as a Python binding for the C
library of Jim Kent's utilties for the UCSC genome browser.


See also
--------
`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`Source repository for Kent utilities <https://github.com/ENCODE-DCC/kentUtils.git>`_
    The header files are particularly useful.
"""
from plastid.readers.bbifile cimport lm, _BBI_Reader




cdef class BigWigReader(_BBI_Reader):
    cdef:
        lm* _lm
        double fill
        
    cdef lm* _get_lm(self)
#    cdef np.ndarray segment_counts(self,GenomicSegment)
#    cdef np.ndarray chain_counts(self,SegmentChain)
    

#cdef class BigBedFile(BBI_File)
#    cdef list what_overlaps_segment(self,GenomicSegment)
#    cdef list what_overlaps_chain(self,SegmentChain)
