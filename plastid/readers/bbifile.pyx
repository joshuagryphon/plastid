"""Base class for reading BBI (BigWig & BigBeD) files. This is a Cython wrapper
for Jim Kent's utilities.
"""
import os

cdef class _BBI_Reader:
    """Abstract base class for BigBed and BigWig file readers

    Reads basic file properties
    """
    def __cinit__(self, str filename, *args, **kwargs):
        if not os.path.exists(filename):
            raise IOError("File '%s' not found." % filename)
        
        self.filename = filename
        self._chrominfo = None
        self._zoomlevels = None
        self.offsets = None

    def __dealloc__(self):
        """Close BBI file"""
        close_file(self._bbifile)

    property kind:
        def __get__(self):
            """Type signature for file"""
            return self._bbifile.typeSig

    property version:
        """Version of BigWig or BigBed file format"""
        def __get__(self):
            return self._bbifile.version

    property chroms:
        """Dictionary mapping chromosome names to lengths"""
        def __get__(self):
            return self.c_chroms()
    
#    property zoomlevels:
#        pass
#    
#    property autosql:
#        pass

    property uncompress_buf_size:
        """Size of buffer needed to uncompress blocks. If 0, the data is uncompressed"""
        def __get__(self):
            return self._bbifile.uncompressBufSize

    cdef dict _define_chroms(self):
        """Return dictionary mapping chromosome names to their lengths

        Returns
        -------
        dict
            Dictionary mapping chromosome names to their lengths
        """
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
        """Return dictionary mapping chromosome names to their lengths,
        populating ``self._chrominfo`` if necessary
        
        Returns
        -------
        dict
            Dictionary mapping chromosome names to their lengths
        """
        if self._chrominfo is not None:
            return self._chrominfo
        else: 
            self._chrominfo = self._define_chroms()
            return self._chrominfo
