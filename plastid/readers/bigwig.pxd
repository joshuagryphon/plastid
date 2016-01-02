from plastid.readers.bbifile cimport bbiFile, bbiChromInfo, \
                                 bbiChromList, bbiChromInfoFreeList, \
                                 bbiCachedChromLookup, bigWigValsOnChrom, \
                                 bigWigFileOpen, bigWigIntervalQuery, close_file,\
                                 _BBI_File


cdef class _BBI_File:

    cdef:
        bbiFile * _bbifile
        str filename
        dict _chrominfo
        dict _zoomlevels
        dict offsets

    cdef dict _define_chroms(self)
    cdef dict c_chroms(self)


cdef class BigWigFile(_BBI_File):
    pass
#    cdef np.ndarray segment_counts(self,GenomicSegment)
#    cdef np.ndarray chain_counts(self,SegmentChain)
    

#cdef class BigBedFile(BBI_File)
#    cdef list what_overlaps_segment(self,GenomicSegment)
#    cdef list what_overlaps_chain(self,SegmentChain)
