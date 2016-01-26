from collections import OrderedDict

from plastid.readers.bbifile cimport bbiFile, bits32, bits64, lm, lmInit, lmCleanup, freeMem
from plastid.readers.bbifile cimport _BBI_Reader
from plastid.util.services.mini2to3 import safe_bytes, safe_str
from plastid.genomics.roitools cimport GenomicSegment, SegmentChain
from plastid.genomics.c_common cimport strand_to_str, str_to_strand, Strand, \
                                       forward_strand, reverse_strand, unstranded,\
                                       error_strand,\
                                       _GeneratorWrapper

cdef extern from "<bigBed.h>":

    cdef struct bigBedInterval:
        bigBedInterval * next
        bits32           start, end
        char           * rest                #remainder of BED line
        bits32           chromId

    bbiFile* bigBedFileOpen(char* fileName)
    
# /* Get data for interval.  Return list allocated out of lm.  Set maxItems to maximum
#  * number of items to return, or to 0 for all items. */
    bigBedInterval * bigBedIntervalQuery(bbiFile *bbi, char *chrom, bits32 start,
                                         bits32 end, int maxItems, lm *lm)
    
    
# /* Convert bigBedInterval into an array of chars equivalent to what you'd get by
#  * parsing the bed file. The startBuf and endBuf are used to hold the ascii representation of
#  * start and end.  Note that the interval->rest string will have zeroes inserted as a side effect.
#  * Returns number of fields in row.  */    
    int bigBedIntervalToRow(bigBedInterval *interval, char *chrom, char *startBuf,
                            char *endBuf, char **row, int rowSize)
    

    # number of items in file
    bits64 bigBedItemCount(bbiFile *bbi)
    
    # get autoSql
    char *bigBedAutoSqlText(bbiFile *bbi)
    
# /* Return index associated with fieldName.  Aborts if no such index.  Optionally return
#  * index in a row of this field. */
#    bptFile *bigBedOpenExtraIndex(bbiFile *bbi, char *fieldName, int*retFieldIX)
    
# /* Return list of names of extra indexes beyond primary chorm:start-end one */
#    slName *bigBedListExtraIndexes(bbiFile *bbi)
 
cdef class BigBedReader(_BBI_Reader):
    cdef:
        bint               add_three_for_stop
        int                total_fields
        int                bed_fields
        int                num_extension_fields
        object             extension_fields  # maps name to description
        object             extension_types   # map  name to formatting func
        
        object    return_type
        #types.classTypes return_type

    cdef _GeneratorWrapper _c_get(self, SegmentChain roi, bint stranded=*, bint check_unique=*)
    cdef str process_record(self,dict chromids,bigBedInterval *iv)
#    cpdef list _getiter(self,GenomicSegment,bint)    