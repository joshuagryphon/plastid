from collections import OrderedDict # cimport?

from libc.stddef cimport size_t
from plastid.genomics.c_common cimport Strand, forward_strand, reverse_strand, unstranded
from plastid.genomics.roitools cimport GenomicSegment

cdef:
    str WARN_CHROM_NOT_FOUND

#===============================================================================
# INDEX: Common info for all bigwig files 
#===============================================================================

DEF bigBedSig = 0x8789F2EB
DEF bigWigSig = 0x888FFC26 

cdef extern from "<common.h>" nogil:
    ctypedef unsigned char   Bits
    ctypedef unsigned char   bits8
    ctypedef unsigned short  bits16
    ctypedef unsigned bits32
    ctypedef unsigned long long bits64

    void freeMem(void *pt)
#     
#     cdef struct fileOffsetSize:
#         fileOffsetSize * next
#         bits64           offset
#         bits64           size
    


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
        bint   isSwapped
        #struct bptFile *chromBpt
        bits16 version
        bits16 zoomLevels
        #bits64 chromTreeOffset
#         bits64 unzoomedDataOffset
#         bits64 unzoomedIndexOffset
        bits16 fieldCount
        bits16 definedFieldCount
        bits64 asOffset
        bits64 totalSummaryOffset
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


# cdef extern from "<cirTree.h>":
#     cdef struct cirTreeFile:
#     # R tree index file handle. */
#         struct cirTreeFile *next;    # Next in list of index files if any. */
#         char *fileName;              # Name of file - for error reporting. */
#         #struct udcFile *udc;         # Open file pointer. */
#         boolean isSwapped;           # If TRUE need to byte swap everything. */
#         bits64 rootOffset;           # Offset of root block. */
#         bits32 blockSize;            # Size of block. */
#         bits64 itemCount;            # Number of items indexed. */
#         bits32 startChromIx;         # First chromosome in file. */
#         bits32 startBase;            # Starting base position. */
#         bits32 endChromIx;           # Ending chromosome in file. */
#         bits32 endBase;              # Ending base position. */
#         bits64 fileSize;             # Total size of index file. */
#         bits32 itemsPerSlot;         # Max number of items to put in each index slot at lowest level. */
#     
#         struct fileOffsetSize *cirTreeEnumerateBlocks(struct cirTreeFile *crf)




# Close a BigBed/BigWig file
cdef inline void close_file(bbiFile *myFile):
    bbiFileClose(&myFile)

# Allocate or destroy and reallocate a local memory block of size <= maxmem
cdef lm * get_lm(lm * my_lm=*, int maxmem=*) except NULL


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
        long      _maxmem

#         dict _zoomlevels
#         dict offsets
#         bits64 _unzoomed_data_offset, _unzoomed_index_offset

    cdef dict _define_chroms(self)
    cdef dict c_chroms(self)
    cdef lm* _get_lm(self)
#     cdef dict fetch_summary(self)
