from plastid.readers.bbifile cimport bbiFile, bbiChromInfo,\
                                 bbiChromList, bbiChromInfoFreeList, \
                                 bbiCachedChromLookup, bigWigValsOnChrom, \
                                 bigWigFileOpen, bigWigIntervalQuery, close_file


cdef class BBI_File:
    """Base class for BigBed and BigWig file readers

    Reads basic file properties

    """

    def __cinit__(self,bytes filename):
        self.filename = str(filename)

        # FIXME: lookup
        self._chrominfo = None
        self._zoomlevels = None
        self.offsets = None

    def __dealloc__(self):
        close_file(self._bbifile)

    property kind:
        def __get__(self):
            return self._bbifile.typeSig

    property version:
        def __get__(self):
            return self._bbifile.version

    property chroms:
        def __get__(self):
            return self.c_chroms()
    
#    property zoomlevels:
#        pass
#    
#    property autosql:
#        pass

    property uncompress_buf_size:
        def __get__(self):
            return self._bbifile.uncompressBufSize

    cdef dict _define_chroms(self):
        cdef:
            dict cdict = {}
            bytes chr_name
            int chr_size
            bbiChromInfo *chrom_info = bbiChromList(self._bbifile)

        while chrom_info is not NULL:
            chr_name = chrom_info.name
            chr_size = chrom_info.size
            cdict[chr_name] = chr_size
            chrom_info = chrom_info.next

        bbiChromInfoFreeList(&chrom_info)

        return cdict

    cdef dict c_chroms(self):
        if self._chrominfo is not None:
            return self._chrominfo
        else:
            self._chrominfo = self._define_chroms()
            return self._chrominfo


cdef class BigWigFile(BBI_File):

    def __cinit__(self,bytes filename):
        self._bbifile = bigWigFileOpen(filename)
  

#cdef class BigBedFile(BBI_File):
#    pass
