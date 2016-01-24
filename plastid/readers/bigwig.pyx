"""Reader for `BigWig`_ files, built atop `Jim Kent's utilities`_.


See also
--------
`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`Source repository for Kent utilities <https://github.com/ENCODE-DCC/kentUtils.git>`_
    The header files are particularly useful.
"""
import warnings
cimport numpy
import numpy

from plastid.util.services.exceptions import DataWarning
from plastid.genomics.roitools cimport GenomicSegment, SegmentChain
from plastid.genomics.c_common cimport reverse_strand
from plastid.readers.bbifile cimport _BBI_Reader, close_file,\
                                     bbiSummary, bbiSummaryType,\
                                     bbiSumMax,bbiSumMin,bbiSumMean,\
                                     bbiSumCoverage, bbiSumStandardDeviation,\
                                     lmInit, lmCleanup, lmAlloc, lm, \
                                     bitFindClear

from plastid.readers.bbifile cimport WARN_CHROM_NOT_FOUND
from plastid.util.services.mini2to3 import safe_bytes, safe_str

#===============================================================================
# INDEX: BigWig reader
#===============================================================================



cdef class BigWigReader(_BBI_Reader):
    """Reader providing random or sequential access to data stored in `BigWig`_ files.
    """

    def __cinit__(self, str filename, **kwargs): #, double fill = numpy.nan):
        """Open a `BigWig`_ file.
        
        Parameters
        ----------
        filename : str
            Name of `bigwig`_ file
        """
#         """
#         fill : float
#             Value to use when there is no data covering a base (e.g. zero, nan,
#             et c. Default: `numpy.nan`)
#         """
        self._bbifile = bigWigFileOpen(safe_bytes(filename))
        self.fill = 0.0 #fill
        self._sum = numpy.nan

    def __iter__(self):
        return BigWigIterator(self)

    def _get(self, GenomicSegment roi, bint roi_order=True, double fill=numpy.nan): 
        """Retrieve array of counts from a region of interest.
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Region of interest in genome
        
        roi_order : bool, optional
            If `True` (default) return vector of values 5' to 3' 
            relative to vector rather than genome.

        fill : double, optional
            Override fill value to put in positions with no data.
            (Default: Use value of `self.fill`)


        Returns
        -------
        :class:`numpy.ndarray`
            vector of numbers, each position corresponding to a position
            in `roi`, from 5' to 3' relative to `roi`

        
        See also
        --------
        plastid.genomics.roitools.SegmentChain.get_counts
            Fetch a spliced vector of data covering a |SegmentChain|
        """
        cdef:
            long start    = roi.start
            long end      = roi.end
            str  chrom    = roi.chrom
            size_t length = end - start
            double usefill = self.fill if numpy.isnan(fill) else fill
            
            numpy.ndarray counts = numpy.full(length,usefill,dtype=numpy.float)
            double [:] view      = counts
            
            lm* buf = self._get_lm()
            bbiInterval* iv
            long segstart, segend

        # FIXME: add segmentchain here
        
        # return empty vector if chromosome is not in BigWig file
        if chrom not in self.c_chroms():
            warnings.warn(WARN_CHROM_NOT_FOUND % (chrom,self.filename),DataWarning)
            return counts
        
        # populate vector
        iv = bigWigIntervalQuery(self._bbifile,safe_bytes(chrom),start,end,buf)
        while iv is not NULL:
            segstart = iv.start - start
            segend = iv.end - start
            view[segstart:segend] = iv.val
            iv = iv.next
        
        if roi.strand == "-" and roi_order == True:
            counts = counts[::-1]

        return counts
                     
    def __getitem__(self, GenomicSegment roi):
        """Retrieve array of counts at each position in `roi`, in `roi`'s 5' to 3' direction
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Region of interest in genome

        Returns
        -------
        :class:`numpy.ndarray`
            vector of numbers, each position corresponding to a position
            in `roi`, from 5' to 3' relative to `roi`

        
        See also
        --------
        plastid.genomics.roitools.SegmentChain.get_counts
            Fetch a spliced vector of data covering a |SegmentChain|
        """
        return self._get(roi)

    def sum(self):
        """Calculate sum of data in `BigWig`_ file"""
        cdef double mysum = 0.0
        if not numpy.isnan(self._sum):
            return self._sum
        else:
            for chrom in self.chroms:
                # n.b. - we calculate the sum from numpy arrays
                # rather than using the stored sum, which is approximate
                mysum += numpy.nansum(self.get_chromosome(chrom))
                
            self._sum = mysum
            
        return mysum
    
    def get_chromosome(self, str chrom):
        """Retrieve values across an entire chromosome more efficiently than using ``bigwig_reader[chromosome_roi]``
        
        Parameters
        ----------
        chrom : str
            Chromosome name
            
        Returns
        -------
        :class:`numpy.ndarray`
            Numpy array of float values covering entire chromosome `chrom`
        """
        cdef:
            lm* buf = self._get_lm()
            bigWigValsOnChrom* vals
            long length, bufsize
            int i = 0
            numpy.ndarray counts, mask
            bint success
            numpy.double_t [:] valview
            double fill = self.fill
            
        #     cdef struct bigWigValsOnChrom:
        #         bigWigValsOnChrom *next
        #         char   *chrom
        #         long   chromSize
        #         long   bufSize     # size of allocated buffer
        #         double *valBuf     # value for each base on chrom. Zero where no data
        #         Bits   *covBuf     # a bit for each base with data
        vals    = bigWigValsOnChromNew()
        success = bigWigValsOnChromFetchData(vals,safe_bytes(chrom),self._bbifile)
        length  = vals.chromSize

        if chrom not in self.c_chroms():
            warnings.warn(WARN_CHROM_NOT_FOUND % (chrom,self.filename),DataWarning)
            counts = numpy.full(length,fill)
        elif success == False:
            warnings.warn("Could not retrieve data for chrom '%s' from file '%s'." % (chrom,self.filename),DataWarning)
            counts = numpy.full(length,fill)
        else:
            valview = <numpy.double_t[:length]> vals.valBuf
            
            # it would be nice not to copy; but if we borrow the memory
            # allocated by the central stack in Jim Kent's source code,
            # bad things might happen.
            counts = numpy.asarray(valview.copy())
            
            # if fill isn't 0, set no-data values in output to self.fill
            # these are in vals.covBuf
            if self.fill != 0:
                i=-1
                while True:
                    i = bitFindClear(vals.covBuf,i+1,length)
                    if i >= length: # no clear bit left
                        break
                    counts[i] = fill

        bigWigValsOnChromFree(&vals)
        return counts

    def summarize(self, GenomicSegment roi):
        """Summarize `BigWig`_ data over the region of interest.
        
        Note
        ----
        These are quickly calculated approximations. More accurate mean, stdev,
        et c can be calculated using ``bigwig_reader[roi].mean()`` et c
         
        Returns
        -------
        tuple of floats
            `(mean,max,min,fraction of bases covered,stdev)` of data in `BigWig`_ 
            file over the region defined by `roi` 
        """
        mean  = self._summarize(roi,bbiSumMean)
        max_  = self._summarize(roi,bbiSumMax)
        min_  = self._summarize(roi,bbiSumMin)
        cov   = self._summarize(roi,bbiSumCoverage)
        stdev = self._summarize(roi,bbiSumStandardDeviation)
         
        return (mean,max_,min_,cov,stdev)
        
    cdef double _summarize(self, GenomicSegment roi, bbiSummaryType type_):
        """Summarize `BigWig`_ data over ROI for a single statistic
         
        Parameters
        ----------
        roi : |GenomicSegment|
             
        type: int in [0...4]
            Type of data to calculate: mean (type_ = 0), max (1), min (2),
            fraction of bases covered (3), stdev (4)
            
        Returns
        -------
        float
            Summarized value
        """
        cdef:
            double retval
            str chrom = roi.chrom
             
        retval = bigWigSingleSummary(self._bbifile,
                                     safe_bytes(chrom),
                                     roi.start,
                                     roi.end,
                                     type_,
                                     self.fill)
         
        return retval

    property fill:
        def __set__(self,double val):
            self.fill = val
            
        def __get__(self):
            return self.fill


#===============================================================================
# INDEX: BigWig iterator
#===============================================================================

####################################################################
# for some reason, the return value for this gives no repr in python
# so:
#
#     >>> z = iter(some_bigwig_file) gives
#     >>> z
#
# gives:
#
# TypeError: object of type 'getset_descriptor' has no len()
#
# while:
#
#    >>>> repr(z)
#
# gives:
#
# '<generator object at 0x7f066e0f3d98>'
#

def BigWigIterator(_BBI_Reader reader):
    """Iterate over values in the `BigWig`_ file, sorted lexically by chromosome and position.
     
    Yields
    ------
    tuple
        `(chrom name, start, end, value)`, where start & end are zero-indexed
        and half-open
        
    Raises
    ------
    MemoryError
        If memory cannot be allocated
    """
    cdef:
        dict chromsizes = reader.c_chroms()
        list chroms = sorted(chromsizes.keys())
        str chrom
        long chromlength
        double val
        lm* buf
        bbiInterval* iv
     
    buf = lmInit(0)
    for chrom in chroms:
        chromlength = chromsizes[chrom]
        iv = bigWigIntervalQuery(reader._bbifile,
                                 safe_bytes(chrom),
                                 0,
                                 chromlength,
                                 buf)
        while iv is not NULL:
            retval = (chrom,long(iv.start),long(iv.end),float(iv.val))
            iv = iv.next
            yield retval

    lmCleanup(&buf)
