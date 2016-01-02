from plastid.genomics.roitools cimport GenomicSegment, SegmentChain
from plastid.genomics.c_common cimport reverse_strand

from plastid.readers.bbifile cimport bbiFile, bbiInterval, bbiChromInfo,\
                                 bbiChromList, bbiChromInfoFreeList, \
                                 bbiCachedChromLookup, bigWigValsOnChrom, \
                                 bigWigFileOpen, bigWigIntervalQuery, close_file,\
                                 lmInit, lmCleanup, lmAlloc, lm

cimport numpy
import numpy


cdef class _BBI_File:
    """Abstract base class for BigBed and BigWig file readers

    Reads basic file properties
    """
    def __cinit__(self, str filename):
        self.filename = filename
        self._chrominfo = None
        self._zoomlevels = None
        self.offsets = None

    def __dealloc__(self):
        """Close BBI file"""
        close_file(self._bbifile)

    property kind:
        def __get__(self):
            """Type signature for file"""
            return self._bbifile.typeSig

    property version:
        """Version of BigWig or BigBed file format"""
        def __get__(self):
            return self._bbifile.version

    property chroms:
        """Dictionary mapping chromosome names to lengths"""
        def __get__(self):
            return self.c_chroms()
    
#    property zoomlevels:
#        pass
#    
#    property autosql:
#        pass

    property uncompress_buf_size:
        """Size of buffer needed to uncompress blocks. If 0, the data is uncompressed"""
        def __get__(self):
            return self._bbifile.uncompressBufSize

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




cdef class BigWigFile(_BBI_File):

    def __cinit__(self,str filename):
        """Open a `BigWig`_ file.
        
        Parameters
        ----------
        filename : str
            Name of bigwig file
        """
        self._bbifile = bigWigFileOpen(bytes(filename))
        # local memory buffer
        self._lm = NULL
        
    def __dealloc__(self):
        """Deallocate memory buffer ``self._lm``, if allocated"""
        if self._lm != NULL:
            lmCleanup(&self._lm)
  
    cdef lm * _get_lm(self):
        """Return ``self._lm``, allocating it if necessary
        
        Returns
        -------
        lm
            local memory buffer
        """
        if self._lm == NULL:
            self._lm = lmInit(0)
            
        return self._lm
         
    def __getitem__(self,roi): #,roi_order=True): 
        """Retrieve array of counts from a region of interest, following
        the mapping rule set by :meth:`~BAMGenomeArray.set_mapping`.
        Values in the vector are ordered 5' to 3' relative to `roi`
        rather than the genome (i.e. are reversed for reverse-strand
        features).
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Region of interest in genome
        
        roi_order : bool, optional
            If `True` (default) return vector of values 5' to 3' 
            relative to vector rather than genome.

        Returns
        -------
        numpy.ndarray
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
            size_t length = end - start
            numpy.ndarray counts = numpy.zeros(length,dtype=numpy.float)
            double [:] view = counts
            lm* buf = self._get_lm()
            bbiInterval* iv
            long segstart, segend
            
        if not buf:
            raise MemoryError("BigWig.__get__: Could not allocate memory.")
        
        # populate vector
        iv = bigWigIntervalQuery(self._bbifile,roi.chrom,start,end,buf)
        while iv is not NULL:
            print("iv: %s\t%s\t%s" % (iv.start,iv.end,iv.val))
            segstart = iv.start - start
            segend = iv.end - start
            view[segstart:segend] = iv.val
            iv = iv.next

#         if roi_order == True and roi.c_strand == reverse_strand:
#             counts = counts[::-1]

        return counts

    def get_chromosome(self,str chrom):
        pass
        
    def __iter__(self):
        """Iterate over values in the `BigWig`_ file, chromosome by chromosome.
        
        Yields
        ------
        tuple
            `(chrom name, start, end, value)`, where start & end are 0-indexed
            and half-open
        """
        pass  

#cdef class BigBedFile(BBI_File):
#    pass
