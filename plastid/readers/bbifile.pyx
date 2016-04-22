"""Base class for `BigBed`_ and `BigWig`_ readers, both of which are implemented
as a Python binding for the C library of Jim Kent's utilties for the UCSC genome
browser.


See also
--------
`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`Source repository for Kent utilities <https://github.com/ENCODE-DCC/kentUtils.git>`_
    The header files are particularly useful.
"""
import os
import warnings

from plastid.readers.autosql import AutoSqlDeclaration
from plastid.util.io.binary import BinaryParserFactory, find_null_bytes
from plastid.util.services.mini2to3 import safe_bytes, safe_str

#===============================================================================
# Exported warnings
#===============================================================================

WARN_CHROM_NOT_FOUND = "No data for chromosome '%s' in file '%s'." 
WARN_FILE_NOT_FOUND = "File '%s' not found."


#===============================================================================
# Functions
#===============================================================================


cdef lm * get_lm(lm * my_lm=NULL, int maxmem=0) except NULL:
    """Allocate or destroy and reallocate a local memory pool, aiming to keep it under `maxmem`
    
    Parameters
    ----------
    my_lm : *lm, optional
        Local memory pool. If NULL, a new pool is allocated. (Default: NULL)
    
    maxmem : int, optional
        Maximum desired size in bytes for memory pool. If the size of `my_lm`
        exceeds `maxmem`, `my_lm` will be freed and a new pool created.
    
    
    Returns
    -------
    lm
        local memory buffer, or NULL if memory cannot be allocated
        
    Raises
    ------
    MemoryError
        If memory cannot be allocated
    """
    cdef:
        int realmaxmem = maxmem >> 10
        
    if my_lm == NULL:
#         print("Initializing LM in get_lm()")
        my_lm = lmInit(0)

    elif maxmem > 0 and lmSize(my_lm) >= maxmem:
#         print("LM Size (%s kb allocated. %s available) exceeded maxmem (%s kb). Cleaning up ..." % (lmSize(my_lm) >> 10,
#                                                                                                     lmAvailable(my_lm) >> 10,
#                                                                                                     realmaxmem))
        
        lmCleanup(&my_lm)
        my_lm = NULL
        my_lm = lmInit(0)
#         print("    New size: %s kb (%s kb available). Max: %s" % (lmSize(my_lm) >> 10,
#                                                                   lmAvailable(my_lm) >> 10,
#                                                                   realmaxmem>>10))

    if not my_lm:
        raise MemoryError("BigBed/BigWig reader: Could not allocate local memory block.")

#     print("LM Size (%s kb allocated. %s available. %s max)." % (lmSize(my_lm) >> 10,lmAvailable(my_lm) >> 10,maxmem>>10))
    return my_lm


#===============================================================================
# Classes
#===============================================================================


cdef class _BBI_Reader:
    """Abstract base class for `BigWig`_ file readers

    Reads basic file properties.

    Parameters
    ----------
    filename : str
        Name of file to open

    maxmem : float
        Maximum desired memory footprint for C objects, in megabytes.
        May be temporarily exceeded if large queries are requested.
        (Default: 0, No maximum)
    """
    def __cinit__(self, str filename, maxmem=0,*args, **kwargs):
        if not os.path.exists(filename):
            raise IOError(WARN_FILE_NOT_FOUND % filename)
        
        self.filename      = filename
        self._chromlengths = None
        self._chromids     = None
        self._summary      = None
        self._lm           = NULL
        self._maxmem       = int(round(maxmem * 1024 * 1024))

    def __dealloc__(self):
        """Close `BigBed`_/`BigWig`_ file"""

        # deallocate memory, if allocated
        if self._lm != NULL:
            lmCleanup(&self._lm)

        # close file
        close_file(self._bbifile)

    def __repr__(self):
        return "<%s file='%s'>" % (self.__class__.__name__,self.filename)

    def __str__(self):
        return repr(self)

    cdef lm * _get_lm(self):
        """Return ``self._lm``, allocating it if necessary
        
        Returns
        -------
        lm
            local memory buffer
            
        Raises
        ------
        MemoryError
            If memory cannot be allocated
        """
        self._lm = get_lm(my_lm=self._lm,maxmem=self._maxmem)
        return self._lm

    property filename:
        """Name of BigWig or BigBed file"""
        def __get__(self):
            return self.filename
        
    property version:
        """Version of BigWig or BigBed file format"""
        def __get__(self):
            return self._bbifile.version

    property chroms:
        """Dictionary mapping chromosome names to lengths"""
        def __get__(self):
            return self.c_chroms()

    property chrom_sizes:
        """DEPRECATED: Use `.chroms` instead of `.chrom_sizes`"""
        def __get__(self):
            warnings.warn("`chrom_sizes` is deprecated and will be removed from plastid v0.6.0. Use `chroms` instead",UserWarning)
            return self.c_chroms()

    # FIXME: hide    
    property chromids:
        def __get__(self):
            return self._chromids
# see note below.
#     property summary_info:
#         """Summary information over total BBI file"""
#         def __get__(self):
#             if self._summary is not None:
#                 return self._summary
#             else:
#                 return self.fetch_summary()

    property uncompress_buf_size:
        """Size of buffer needed to uncompress blocks. If 0, the data is uncompressed"""
        def __get__(self):
            return self._bbifile.uncompressBufSize

# data for mean as read from table is very inaccurate. Decide what to do later    
#     cdef dict fetch_summary(self):
#         """Return summary info of BBI file, and set `self._summary`.
#         Summary info includes:
#         
#           - min value
#           - max value
#           - sum of values
#           - sum of squares of values
#           - covered bases
#           - mean value over covered bases
#           
#         Returns
#         -------
#         dict
#             dictionary described above
#         """
#         cdef:
#             bbiSummaryElement mydata = bbiTotalSummary(self._bbifile)
#             dict dtmp
#          
#         dtmp = {
#             "min" : mydata.minVal,
#             "max" : mydata.maxVal,
#             "mean" : mydata.sumData / mydata.validCount,
#             "sum"  : mydata.sumData,
#             "sum_squares" : mydata.sumSquares,
#             "valid_bases" : mydata.validCount,
#         }
#         self._summary = dtmp
#         return dtmp

    cdef dict _define_chroms(self):
        """Return dictionary mapping chromosome names to their lengths

        Returns
        -------
        dict
            Dictionary mapping chromosome names to their lengths
        """
        cdef:
            dict cdict = {}
            dict idict = {}
            str chr_name
            int chr_size
            bbiChromInfo *chrom_info = bbiChromList(self._bbifile)

        while chrom_info is not NULL:
            chr_name         = safe_str(chrom_info.name)
            chr_size         = chrom_info.size
            cdict[chr_name]  = chr_size
            idict[chrom_info.id]        = chr_name
            chrom_info       = chrom_info.next

        self._chromlengths = cdict
        self._chromids     = idict 
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
        if self._chromlengths is not None:
            return self._chromlengths
        else: 
            self._define_chroms()
            return self._chromlengths

