"""Base class for `BigWig`_ reader, implemented as a Python binding for the C
library of Jim Kent's utilties for the UCSC genome browser.


See also
--------
`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`Source repository for Kent utilities <https://github.com/ENCODE-DCC/kentUtils.git>`_
    The header files are particularly useful.
"""

from collections import OrderedDict # cimport?

from libc.stddef cimport size_t
from plastid.genomics.c_common cimport Strand, forward_strand, reverse_strand, unstranded
from plastid.genomics.roitools cimport GenomicSegment

cdef:
    str WARN_CHROM_NOT_FOUND

#===============================================================================
# INDEX: Common info for all bigwig filefs 
#===============================================================================

DEF bigBedSig = 0x8789F2EB
DEF bigWigSig = 0x888FFC26 

cdef extern from "<common.h>" nogil:
    ctypedef unsigned char Bits
    ctypedef unsigned char  bits8
    ctypedef unsigned short bits16
    ctypedef unsigned bits32
    ctypedef unsigned long long bits64

cdef extern from "<bits.h>":
    
    # find index of next set bit
    int bitFindSet(Bits *b, int startIx, int bitCount)
    
    # find index of next clear bit
    int bitFindClear(Bits *b, int startIX, int bitCount)
    
    # read next bit
    bint bitReadOne(Bits *b, int bitIx)

cdef extern from "<localmem.h>":
    cdef struct lm:
        lmBlock *blocks
        size_t blockSize
        size_t allignMask
        size_t allignAdd
        
    cdef struct lmBlock:
        lmBlock *next
        char *free
        char *end
        char *extra
    
    lm * lmInit(int blockSize)
    size_t lmAvailable(lm *lm)
    size_t lmSize(lm *lm)
    void * lmAlloc(lm *lm, size_t size)
    void * lmAllocMoreMem(lm *lm, void *pt, size_t oldSize, size_t newSize)
    void lmCleanup(lm **pLm)

cdef extern from "<bbiFile.h>":
    cdef struct bbiZoomLevel:
        bbiZoomLevel *next 
        bits32 reductionLevel
        bits64 dataOffset
        bits64 indexOffset

    cdef struct bbiFile:
        bbiFile *next
        char *fileName
        #struct udcFile *udc
        bits32 typeSig
        #struct bptFile *chromBpt
        bits16 version
        bits16 zoomLevels
        #bits64 chromTreeOffset
#         bits64 unzoomedDataOffset
#         bits64 unzoomedIndexOffset
        bits16 fieldCount
        bits16 definedFieldCount
        bits64 asOffset
        #bits64 totalSummaryOffset
        bits32 uncompressBufSize
        #bits64 extensionOffset
        #struct cirTreeFile unzoomedCir
        bbiZoomLevel *levelList
        #bits16 extensionSize
        bits16 extraIndexCount
        #bits64 extraIndexListOffset

    cdef struct bbiChromIdSize:
        bits32 chromId
        bits32 chromSize

    cdef struct bbiChromInfo:
        bbiChromInfo * next
        char *name
        bits32 id
        bits32 size

    cdef struct bbiSummary:
        bbiSummary * next_
        bits32 chromId
        bits32 start, end
        bits32 validCount
        float minVal
        float maxVal
        float sumData
        float sumSquares
        bits64 fileOffset

    cdef struct bbiSummaryElement:
        bits64 validCount
        double minVal
        double maxVal
        double sumData
        double sumSquares

    cdef enum bbiSummaryType:
        bbiSumMean              = 0
        bbiSumMax               = 1
        bbiSumMin               = 2
        bbiSumCoverage          = 3
        bbiSumStandardDeviation = 4

    # bbiFile *bbiFileOpen(char *fileName, bits32 sig, char *typeName)

    void bbiFileClose(bbiFile **pBwf)

    bbiChromInfo *bbiChromList(bbiFile *bbi)
    void bbiChromInfoFreeList(bbiChromInfo **pList)

    char *bbiCachedChromLookup(bbiFile *bbi,
                               int chromId,
                               int lastChromId,
                               char *chromBuf,
                               int chromBufSize)

    bbiSummaryElement bbiTotalSummary(bbiFile *bbi)



cdef inline void close_file(bbiFile *myFile):
    bbiFileClose(&myFile)



# cdef extern from "bigBed.h":
# 
#     cdef struct bigBedInterval:
#         bigBedInterval *next
#         bits32 start, end
#         char *rest
#         bits32 chromId
# 
#     bbiFile *bigBedFileOpen(char *fileName)
#     bigBedInterval bigBedIntervalQuery(bbiFile *bbi, char *chrom, bits32 start, bits32 end,
#                                        int maxItems, struct lm *lm)
# 
#     bint bigBedSummaryArray(bbiFile *bbi, char *chrom, bits32 start, bits32 end,
#                             enum bbiSummaryType summaryType, int summarySize, double *summaryValues)
# 
# #    int bigBedIntervalToRow(bigBedInterval *interval, char *chrom, char *startBuf, char*endBuf,
# #                            char **row, int rowSize)
# #
# #    int bigBedIntervalToRowLookupChrom(bigBedInterval *interval,
# #                                       bigBedInterval *prevInterval,
# #                                       bbiFile *bbi,
# #                                       char *chromBuf,
# #                                       int chromBufSize,
# #                                       char *startBuf,
# #                                       char *endBuf,
# #                                       char **row,
# #                                       int rowSize)
# #
# ##    struct bigBedInterval *bigBedNameQuery(bbiFile *bbi, bptFile *index, int fieldIx, char *name, lm *lm)
# #
# #    struct bigBedInterval *bigBedMultiNameQuery(bbiFile *bbi, bptFile *index, int fieldIx, char **names,
# #                                                int nameCount, lm * lm)'
# #
#     bits64 bigBedItemCount(bbiFile *bbi)
#     char * bigBedAutoSqlText(bbiFile *bbi)


#===============================================================================
# INDEX: Python class
#===============================================================================


cdef class _BBI_Reader:

    cdef:
        bbiFile * _bbifile
        str       filename
        dict      _chromids
        dict      _chromlengths
        dict      _summary
        lm *      _lm

#         dict _zoomlevels
#         dict offsets
#         bits64 _unzoomed_data_offset, _unzoomed_index_offset

    cdef dict _define_chroms(self)
    cdef dict c_chroms(self)
    cdef lm* _get_lm(self)
#     cdef dict fetch_summary(self)


# cdef class bb2(_BBI_Reader):
#     cdef:
#         _lm          * lm
#         class          return_type
#         OrderedDict    custom_fields
#         bits64         num_records
#         bint           add_three_for_stop
# 
#     cdef list c_getitem(self,GenomicSegment)
#     cdef str c_get_autosql(self)
#     cdef *lm _get_lm(self)

