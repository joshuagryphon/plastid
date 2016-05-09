"""Fast reader for `BigWig`_ files, implemented as a Python binding for the C
library of Jim Kent's utilties for the UCSC genome browser.


See also
--------
`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`Source repository for Kent utilities <https://github.com/ENCODE-DCC/kentUtils.git>`_
    The header files are particularly useful.
"""
from plastid.readers.bbifile cimport lm, bbiFile, _BBI_Reader, bits32, Bits, bbiSummaryType, get_lm
from plastid.genomics.roitools cimport GenomicSegment
from plastid.genomics.c_common cimport _GeneratorWrapper

#===============================================================================
# Externs from Kent utilties
#===============================================================================

cdef extern from "<bigWig.h>":
    cdef struct bigWigValsOnChrom:
        bigWigValsOnChrom *next
        char   *chrom
        long   chromSize
        long   bufSize     # size of allocated buffer
        double *valBuf     # value for each base on chrom. Zero where no data
        Bits   *covBuf     # a bit for each base with data

    cdef struct bbiInterval:
        bbiInterval * next
        bits32 start, end
        double val

    bbiFile     * bigWigFileOpen(char *fileName)
    
    bbiInterval * bigWigIntervalQuery(bbiFile *bwf,
                                      char *chrom,
                                      bits32 start,
                                      bits32 end,
                                      lm *lm)
    
    bigWigValsOnChrom *bigWigValsOnChromNew()
    
    void bigWigValsOnChromFree(bigWigValsOnChrom **pChromVals)
    
    bint bigWigValsOnChromFetchData(bigWigValsOnChrom *chromVals,
                                    char *chrom,
                                    bbiFile *bigWig)

    double bigWigSingleSummary(bbiFile *bwf, char *chrom, int start, int end,
                               bbiSummaryType summaryType, double defaultVal)


#===============================================================================
# Python class
#===============================================================================


cdef class BigWigReader(_BBI_Reader):
    cdef double fill
    cdef double _sum
    
    cdef double _summarize(self,GenomicSegment roi, bbiSummaryType type_)
    cdef double c_sum(self)
    cdef bigWigValsOnChrom * c_get_chromosome_counts(self, str chrom)

