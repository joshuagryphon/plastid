"""Base class for `BigWig`_ reader, implemented as a Python binding for the C
library of Jim Kent's utilties for the UCSC genome browser.


See also
--------
`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`Source repository for Kent utilities <https://github.com/ENCODE-DCC/kentUtils.git>`_
    The header files are particularly useful.
"""

from libc.stddef cimport size_t


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
        #bits64 exgtensionOffset
        #struct cirTreeFile unzoomedCir
        #bbiZoomLevel *levelList
        #bits16 extensionSize
        #bits16 extraIndexCount
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


#===============================================================================
# INDEX: Python class
#===============================================================================


cdef class _BBI_Reader:

    cdef:
        bbiFile * _bbifile
        str filename
        dict _chrominfo
        dict _summary
#         dict _zoomlevels
#         dict offsets
#         bits64 _unzoomed_data_offset, _unzoomed_index_offset

    cdef dict _define_chroms(self)
    cdef dict c_chroms(self)
#     cdef dict fetch_summary(self)
