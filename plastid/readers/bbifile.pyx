"""Base class for `BigWig`_ reader, implemented as a Python binding for the C
library of Jim Kent's utilties for the UCSC genome browser.


See also
--------
`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`Source repository for Kent utilities <https://github.com/ENCODE-DCC/kentUtils.git>`_
    The header files are particularly useful.
"""
import os


#===============================================================================
# Exported warnings
#===============================================================================

WARN_CHROM_NOT_FOUND = "No data for chromosome '%s' in file '%s'." 
WARN_FILE_NOT_FOUND = "File '%s' not found."


#===============================================================================
# Classes
#===============================================================================

cdef class _BBI_Reader:
    """Abstract base class for BigBed and BigWig file readers

    Reads basic file properties
    """
    def __cinit__(self, str filename, *args, **kwargs):
        if not os.path.exists(filename):
            raise IOError(WARN_FILE_NOT_FOUND % filename)
        
        self.filename = filename
        self._chrominfo = None
        self._summary = None

    def __dealloc__(self):
        """Close BBI file"""
        close_file(self._bbifile)

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
