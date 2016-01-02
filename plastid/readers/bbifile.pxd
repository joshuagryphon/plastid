from libc.stddef cimport size_t
from mercurial.bdiff import blocks
from _sha import blocksize

#cimport numpy as np
#from plastid.genomics.roitools cimport GenomicSegment, SegmentChain

DEF bigBedSig = 0x8789F2EB
DEF bigWigSig = 0x888FFC26 

cdef extern from "<common.h>" nogil:
    ctypedef unsigned char Bits
    ctypedef unsigned char  bits8
    ctypedef unsigned short bits16
    ctypedef unsigned bits32
    ctypedef unsigned long long bits64


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
        #bits64 unzoomedDataOffset
        #bits64 unzoomedIndexOffset
        bits16 fieldCount
        bits16 definedFieldCount
        bits64 asOffset
        #bits64 totalSummaryOffset
        bits32 uncompressBufSize
        #bits64 exgtensionOffset
        #struct cirTreeFile unzoomedCir
        bbiZoomLevel *levelList

        bits16 extensionSize
        bits16 extraIndexCount
        bits64 extraIndexListOffset

    cdef struct bbiChromIdSize:
        bits32 chromId
        bits32 chromSize

    cdef struct bbiChromInfo:
        bbiChromInfo * next
        char *name
        bits32 id
        bits32 size

#    cdef struct bbiSummary:
#        bbiSummary * next_
#        bits32 chromId
#        bits32 start, end
#        bits32 validCount
#        float minVal
#        float maxVal
#        float sumData
#        float sumSquares
#        bits64 fileOffset
#
    cdef struct bbiInterval:
        bbiInterval * next
        bits32 start, end
        double val

    # bbiFile *bbiFileOpen(char *fileName, bits32 sig, char *typeName)

    void bbiFileClose(bbiFile **pBwf)

    bbiChromInfo *bbiChromList(bbiFile *bbi)
    void bbiChromInfoFreeList(bbiChromInfo **pList)

    char *bbiCachedChromLookup(bbiFile *bbi,
                               int chromId,
                               int lastChromId,
                               char *chromBuf,
                               int chromBufSize)

cdef extern from "<bigWig.h>":
    bbiFile     * bigWigFileOpen(char *fileName)
    bbiInterval * bigWigIntervalQuery(bbiFile *bwf,
                                      char *chrom,
                                      bits32 start,
                                      bits32 end,
                                      lm *lm)
    #boolean isBigWig(char *fileName)

    cdef struct bigWigValsOnChrom:
        bigWigValsOnChrom *next
        char   *chrom
        long   chromSize
        long   bufSize     # size of allocated buffer
        double *valBuf     # value for each base on chrom. Zero where no data
        Bits   *covBuf     # a bit for each base with data


cdef inline void close_file(bbiFile *myFile):
    bbiFileClose(&myFile)


