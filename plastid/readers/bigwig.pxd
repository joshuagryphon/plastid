from plastid.readers.bbifile cimport lm, _BBI_File




cdef class BigWigFile(_BBI_File):
    cdef:
        lm* _lm
        
    cdef lm* _get_lm(self)
#    cdef np.ndarray segment_counts(self,GenomicSegment)
#    cdef np.ndarray chain_counts(self,SegmentChain)
    

#cdef class BigBedFile(BBI_File)
#    cdef list what_overlaps_segment(self,GenomicSegment)
#    cdef list what_overlaps_chain(self,SegmentChain)
