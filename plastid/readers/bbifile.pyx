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

from plastid.readers.autosql import AutoSqlDeclaration
from plastid.util.io.binary import BinaryParserFactory, find_null_bytes

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
        self._lm = NULL

    def __dealloc__(self):
        """Close BBI file"""

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
        if self._lm == NULL:
            self._lm = lmInit(0)

        if not self._lm:
            raise MemoryError("BigWig.__get__: Could not allocate memory.")
            
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



# from bigBed.h
# To stream through whole file
 #    struct bbiFile *bbi = bigBedFileOpen(fileName);
 #    struct bbiChromInfo *chrom, *chromList = bbiChromList(bbi);
 #    for (chrom = chromList; chrom != NULL; chrom = chrom->next)
 #        {
 #        struct lm *lm = lmInit(0);
 #        struct bigBedInterval *list = bigBedIntervalQuery(bbi,chrom->name,0,chrom->size,0,lm);
 #        struct bigBedInterval *el;
 #        for (el = list; el != NULL; el = el->next)
 #            // do something involving chrom, el->start, el->end
 #        lmCleanup(&lm);
 #        }
 #    bigBedFileClose(&bbi);


# def BigBedIterator(bb2 bigbed):
#     """Iterate over `BigBed`_ files by chromosome in lexical order
# 
#     Parameters
#     ----------
#     bb2 : bb2
#         `BigBed` reader
# 
#     Yields
#     ------
#     object
#         object of type `bb2.return_type`
#     """
#     cdef:
#         str    chrom
#         bits32 length
#         list   ltmp
# 
#     for chrom, length in sorted(self.c_chroms().items()):
#         ltmp = bb2[GenomicSegment(chrom,length,0,".")]
#         for item in ltmp:
#             yield item
# 
# 
# cdef class bb2(_BBI_Reader)
# 
#     # FIXME: should be object type?
#     def __cinit__(self, str filename, class return_type = None, str add_three_for_stop = False):
#         cdef str asql
# 
#         self._bbifile = bigBedFileOpen(bytes(filename))
#         self.return_type = return_type
#         self.custom_fields = OrderedDict()
#         self.num_records = bigBedItemCount(self._bbifile)
# 
#         self.return_type = SegmentChain if return_type is None else return_type
#         self.add_three_for_stop = add_three_for_stop
# 
#         # do something with custom fields here
#         if self._bbifile.extraIndexCount > 0:
#             asql = bigBedAutoSqlText(self._bbifile)
#             try:
#                 self.autosql_parser = AutoSqlDeclaration(asql)
#                 self.custom_fields  = OrderedDict(list(self.autosql_parser.field_comments.items())[-self._num_custom_fields:])
#             except AttributeError:
#                 warnings.warn("Could not find or could not parse autoSql declaration in BigBed file '%s': %s" % (self.filename,autosql),FileFormatWarning)
#                 self.custom_fields  = OrderedDict([("custom_%s" % X,"no description") for X in range(self._num_custom_fields)])
#                 self.autosql_parser = lambda x: OrderedDict(zip(self.custom_fields,x.split("\t")[-self._num_custom_fields:]))
# 
#     def __iter__(self):
#         return BigBedIterator(self)
# 
#     def __getitem__(self,GenomicSegment roi):
#         return self.c_getitem(roi)
# 
#     cdef list c_getitem(self, GenomicSegment roi, bint stranded = True):
#         """Return a list of strings from BigBed file, covering `roi`
#         """
#         cdef:
#             bigBedInterval iv
#             Strand strand = iv.c_strand
#             list ltmp
#             string stmp
#             list lout = []
#             lm* = self._get_lm()
#             str chrom, stmp, rest
# 
#         # FIXME: deal with strand info
#         iv = bigBedIntervalQuery(self._bbifile,roi.start,roi.end,0,lm)
#         while iv != NULL:
# 
#             # FIXME - get chrom name from ID
#             chrom = "" # some_func(iv.chromId)
# 
#             stmp = "%s\t%s\t%s\t" % (chrom,str(iv.start),str(iv.end)) + iv.rest
# 
#             # TODO: process return type & extra attributes
#             lout.append(stmp)
#             iv = iv.next
# 
#         #FIXME : free iv?
# 
#         return lout
# 
#     cdef str c_get_autosql(self):
#         return bigBedAutoSqlText(self._bbifile)
