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
