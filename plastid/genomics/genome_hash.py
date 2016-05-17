"""This module contains tools for lookup of features in a region of interest within the a genome.

.. contents::
   :local:

Summary
-------

It is frequently useful to retrieve features that overlap specific regions 
of interest in the genome. ``GenomeHashes`` index features by location,
providing quick lookup.


Module contents
---------------

Several implementations are provided, depending how the data are formatted:

====================    =========================================================
**Implementation**      **Format of feature data**
--------------------    ---------------------------------------------------------
|GenomeHash|            Objects in memory or in unindexed `BED`_, `GTF2`_,
                        `GFF3`_, or `PSL`_ files

|BigBedGenomeHash|      Annotations in `BigBed`_ files

|TabixGenomeHash|       Annotations in `tabix`_-compressed `BED`_, `GTF2`_,
                        `GFF3`_, or `PSL`_ files
====================    =========================================================


Examples
--------
Create a |GenomeHash|::

    >>> from plastid import *

    # from objects in memory
    >>> one_hash = GenomeHash(list_of_transcripts)
    
    # from a non-indexed file
    >>> my_hash = GenomeHash(list(GFF3_Reader("some_file.gff")))
    
    # from a BigBed file
    >>> bigbed_hash = BigBedGenomeHash("some_file.bb")
    
    # from tabix-compressed BED file
    >>> tabix_hash = TabixGenomeHash("some_file.bed.gz","BED")


To find features overlapping a region of interest, pass the feature coordinates
to a |GenomeHash| as a |GenomicSegment|, |SegmentChain|, or |Transcript|::

    >>> overlapping = my_hash[GenomicSegment("chrII",50,10000,"+")]
    >>> overlapping
    [ list of SegmentChains / Transcripts, et c ]
    
    # SegmentChains & Transcripts can also be keys:
    >>> tx = Transcript(GenomicSegment("chrII",50,300,"+"),GenomicSegment("chrII",9000,10000,"+"))
    >>> overlapping2 = my_hash[tx]
    >>> overlapping2
    [ list of SegmentChains / Transcripts, et c ]
    
    # find features that overlap `roi` on either strand
    >>> either_strand_overlap = my_hash.get_overlapping_features(roi,stranded=False)
    
"""

#===============================================================================
# Memory-efficient ways to hash features across a genome
#===============================================================================
import copy
from plastid.util.services.mini2to3 import cStringIO
from plastid.util.io.openers import NullWriter, multiopen
from plastid.readers.bed import BED_Reader
from plastid.readers.gff import GTF2_Reader, GFF3_Reader
from plastid.readers.psl import PSL_Reader
from plastid.genomics.roitools import GenomicSegment, SegmentChain
from abc import abstractmethod

DEFAULT_BIN_SIZE=20000


class AbstractGenomeHash(object):
    """Abstract base class for objects that index |SegmentChains| by genomic position,
    allowing quick lookup for comparisons of overlap or other behavior
    """
    @abstractmethod
    def get_overlapping_features(self,roi,stranded=True):
        """Return list of features that overlap `roi`
        
        Parameters
        ----------
        feature : |GenomicSegment| or |SegmentChain|
            Query feature

        stranded : bool
            if `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands (Default: `True`)
                             
                             
        Returns
        -------
        list
            List of overlapping features.


        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        pass

    @abstractmethod
    def __getitem__(self,roi):
        """Return list of features that overlap a region of interest (`roi`)
        on the same strand
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

                             
        Returns
        -------
        list
            Overlapping features


        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        pass


class GenomeHash(AbstractGenomeHash):
    """Index memory-resident features (e.g. |SegmentChains| or |Transcripts|) by genomic position for quick lookup later.
        
     Parameters
     ----------
     features : dict or list, optional
         dict or list of features, as |SegmentChain| objects or subclasses
         (Default: `[]`)
        
     binsize : int, optional
         Size in nucleotides of *neighborhood* for hash. (Default: %s)
        
     do_copy : bool
         If `True`, features will be copied before being stored
         in the hash. This comes at a speed cost, but will prevent unexpected
         side effects if the features are being changed outside the hash. 
         
         If `False` (default), creation of the |GenomeHash| will be much faster.   
         
          
    Notes
    -----
    Because all features are stored in memory, for large genomes, a |TabixGenomeHash|
    or |BigBedGenomeHash| is much more memory-efficient. 
    """
    def __init__(self,features=None,binsize=DEFAULT_BIN_SIZE,do_copy=False):
        """Create a |GenomeHash|
        
         Parameters
         ----------
         features : dict or list, optional
             dict or list of features, as |SegmentChain| objects or subclasses
             (Default: `[]`)
            
         binsize : int, optional
             Size in nucleotides of *neighborhood* for hash. (Default: %s)
            
         do_copy : bool
             If `True`, features will be copied before being stored
             in the hash. This comes at a speed cost, but will prevent unexpected
             side effects if the features are being changed outside the hash. 
             
             If `False` (default), creation of the |GenomeHash| will be much faster.
        """ % DEFAULT_BIN_SIZE
        self.copy = do_copy
        self.feature_dict = {}
        self._id_to_names = {}
        self.binsize = binsize
        features = [] if features is None else features
        self.update(features)
    
    def __repr__(self):
        return "<%s features=%s binsize=%s chrs=%s>" % (self.__class__.__name__,
                                                        len(self.feature_dict),
                                                        self.binsize,
                                                        self._feature_hash.keys())

    def __str__(self):
        return repr(self)
    
    def update(self,features):
        """Add features to the |GenomeHash|
        
        Parameters
        ----------
         features : dict or list
             dict or list of features, as |SegmentChain| objects or subclasses
        """
        try:
            start = 1 + max(self._id_to_names.keys())
        except ValueError:
            # if empty dictionary
            start = 0
        
        if isinstance(features,dict):
            for n,(feature_name,feature) in enumerate(features.items()):
                feature_id = start + n
                self._id_to_names[feature_id] = feature_name
                if self.copy == True:
                    self.feature_dict[feature_id] = copy.deepcopy(feature)
                else:
                    self.feature_dict[feature_id] = feature
        else:
            for n,feature in enumerate(features):
                feature_id = start + n
                feature_name = feature.get_name()
                self._id_to_names[feature_id] = feature_name
                if self.copy == True:
                    self.feature_dict[feature_id] = copy.deepcopy(feature)
                else:
                    self.feature_dict[feature_id] = feature
        
        self._feature_hash = self._make_hash()
        
    def _make_hash(self):
        """Create the |GenomeHash|
           
        Returns
        -------
        dict
            Hierarchichal dictionary: `dict[chrom][strand] = list<str>`
        """
        my_hash = {}
        for feature_id, feature in self.feature_dict.items():
            bins   = self._get_hash_bins(feature)
            chrom  = feature.spanning_segment.chrom
            strand = feature.spanning_segment.strand
            if chrom not in my_hash.keys():
                my_hash[chrom] = {}
                my_hash[chrom]["+"] = {}
                my_hash[chrom]["-"] = {}
            for b in bins:
                try:
                    my_hash[chrom][strand][b].append(feature_id)
                except KeyError:
                    my_hash[chrom][strand][b] = [feature_id]
        return my_hash
        
    def _get_hash_bins(self,roi):
        """Returns a list of genome bins that a given roi falls into
        
        Parameters
        ----------
        feature : |GenomicSegment| or |SegmentChain|
            Query feature
           
           
        Returns
        -------
        list<int>
            List of hash bins that `roi` covers
               
           
        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        if isinstance(roi,GenomicSegment):
            rois = [roi]
        elif isinstance(roi,SegmentChain):
            rois = roi
        else:
            raise TypeError("Query feature must be a GenomicSegment or SegmentChain")
        bins = []
        for seg in rois:
            bin_start = seg.start // self.binsize
            bin_end = seg.end // self.binsize
            bins.extend(range(bin_start,bin_end+1))

        return list(set(bins))
    
    def _get_nearby_feature_ids(self,roi,stranded=True):
        """Return unique IDs of |SegmentChain| s in all the bins occupied by `roi`
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
           Query feature
           
        stranded : bool
            If `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
           
                             
        Returns
        -------
        list<str>
           Unique IDs of features near (within *self.binsize* distance from) `roi`
           
           
        Raises
        ------
        TypeError
            if `roi` is not |GenomicSegment| or |SegmentChain|
        """
        if isinstance(roi,GenomicSegment):
            iv = roi
        elif isinstance(roi,SegmentChain):
            iv = roi.spanning_segment
        else:
            raise TypeError("Query feature must be a GenomicSegment or SegmentChain")
        feature_ranges = self._get_hash_bins(roi)
        nearby_feature_ids = []
        try:
            nearby_hash = self._feature_hash[iv.chrom][iv.strand]
        except KeyError:
            nearby_hash = {}
        
        for b in feature_ranges:
            try:
                nearby_feature_ids.extend(nearby_hash[b])
            except KeyError:
                pass
        
        if stranded is False:
            other_strand = "+" if iv.strand == "-" else "-" if iv.strand == "+" else "." #FIXED
            try:
                nearby_hash = self._feature_hash[iv.chrom][other_strand]
            except KeyError:
                nearby_hash = {}
        
            for b in feature_ranges:
                try:
                    nearby_feature_ids.extend(nearby_hash[b])
                except KeyError:
                    pass
        
        nearby_feature_ids = set(nearby_feature_ids)
        return nearby_feature_ids

    def get_nearby_features(self,roi,stranded=True):
        """Return list of features in all the bins occupied by `roi`
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature

        stranded : bool
            if `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list<str>
           Features near (within `self.binsize` distance from) `roi`
               

        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        nearby_feature_ids = self._get_nearby_feature_ids(roi,stranded=stranded)
        return [self.feature_dict[X] for X in nearby_feature_ids]

    def get_nearby_feature_names(self,roi,stranded=True):
        """Return list of the names of features in all the bins occupied by `roi`
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature

        stranded : bool
            if `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
           Names of features near (within `self.binsize` distance from) `roi`

        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        nearby_feature_ids = self._get_nearby_feature_ids(roi,stranded=stranded)
        return [self._id_to_names[X] for X in nearby_feature_ids]

    def get_overlapping_features(self,roi,stranded=True):
        """Return list of features overlapping `roi`.
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
           Features overlapping `roi`

        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        nearby_features = self.get_nearby_features(roi, stranded=stranded)
        if isinstance(roi,GenomicSegment):
            roi = SegmentChain(roi)

        if stranded == False:
            fn = roi.unstranded_overlaps
        else:
            fn = roi.overlaps
        return [X for X in nearby_features if fn(X) == True]
    
    def __getitem__(self,roi):
        """Return list of features that overlap a region of interest (roi),
        on same strand.
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
           Features overlapping `roi`


        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| |SegmentChain|
        """
        return self.get_overlapping_features(roi,stranded=True)
                

class BigBedGenomeHash(AbstractGenomeHash):
    """BigBedGenomeHash(*filenames,return_type=SegmentChain)
    
    Find features overlapping query regions in `BigBed`_ files.
    
        
    Parameters
    ----------
    *filenames : str 
        One or more filenames to open (NOT open filehandles)

    return_type : class implementing a :py:meth:`from_bed` method
        Class of object to return (Default: |SegmentChain|)
        
    
    Attributes
    ----------
    bigbedreaders : |BigBedReader|
       |BigBedReaders| connecting to BigBed file(s) 
    """
    
    def __init__(self,*filenames,**kwargs): #,base_record_format="III",return_type=None,cache_depth=5):
        """Create a |BigBedGenomeHash|
        
        Parameters
        ----------
        *filenames : str 
            One or more filenames to open (NOT open filehandles)

        return_type : class implementing a :py:meth:`from_bed` method
            Class of object to return (Default: |SegmentChain|)
        """
        from plastid.readers.bigbed import BigBedReader
        return_type = kwargs.get("return_type",SegmentChain)
        
        filenames = list(multiopen(filenames))
        for filename in filenames:
            if not isinstance(filename,str):
                raise ValueError("`filename` must be a 'str'. Found a '%s'." % type(filename)) 
        
        self.filenames = filenames
        self.bigbedreaders = [BigBedReader(X,return_type=return_type) for X in filenames]
    
    def get_overlapping_features(self,roi,stranded=True):
        """Return list of features overlapping `roi`
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            If `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
            Features that overlap `roi`

        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        ltmp = []
        for reader in self.bigbedreaders:
            ltmp.extend(reader.get(roi,stranded=stranded))
            
        return ltmp

    def __getitem__(self,roi,stranded=True):
        """Return list of features that overlap the region of interest (roi)
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            If `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
            Features that overlap `roi`

        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        return self.get_overlapping_features(roi,stranded=stranded)
    

class TabixGenomeHash(AbstractGenomeHash):
    """
    TabixGenomeHash(*filenames,data_format='GTF2')
    
    Find features overlapping query regions in `Tabix`_-indexed files.
        
    Parameters
    ----------
    filenames : str or list of str
        Filename or list of filenames of `Tabix`_-compressed files

    data_format : str
        Format of tabix-compressed file(s). Choices are:
        `'GTF2'`,`'GFF3'`,`'BED'`,`'PSL'` (Default: `GTF2`)
    
    
    Attributes
    ----------
    filenames : str or list
        Name of file to open or list of filenames to open (NOT open filehandles)
        
    tabix_readers : list of :py:class:`pysam.Tabixfile`
       `Pysam`_ interfaces to underlying data files 
    """
    
    _READERS = { "GTF2" : GTF2_Reader,
                 "GFF3" : GFF3_Reader,
                 "BED"  : BED_Reader,
                 "PSL"  : PSL_Reader,
                }
    
    def __init__(self,*filenames,**kwargs): #data_format=None,printer=None):
        """Create a |BigBedGenomeHash|
        
        Parameters
        ----------
        filenames : str or list of str
            Filename or list of filenames of `Tabix`_-compressed files

        data_format : str
            Format of tabix-compressed file(s). Choices are:
            `'GTF2'`,`'GFF3'`,`'BED'`,`'PSL'` (Default: `GTF2`)
        """
        from pysam import Tabixfile
        if len(filenames) == 1 and isinstance(filenames[0],list):
            filenames = filenames[0]
            
        self.filenames = list(multiopen(filenames))
        self.printer = kwargs.get("printer",NullWriter())
        data_format  = kwargs.get("data_format","GTF2")
        try:
            self._reader_class = TabixGenomeHash._READERS[data_format]
        except ValueError:
            msg = "Supported file formats for TabixGenomeHash are: %s" % ", ".join(sorted(TabixGenomeHash._READERS.keys()))
            self.printer.write(msg)
            raise ValueError(msg)
        
        self.tabix_readers = [Tabixfile(X) for X in self.filenames]
    
    def __del__(self):
        try:
            for reader in self.tabix_readers:
                try:
                    reader.close()
                except:
                    pass
        except:
            pass
            
    def get_overlapping_features(self,roi,stranded=True):
        """Return list of features overlapping `roi`
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            If `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
        Returns
        -------
        list
           Features that overlap `roi`

        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        return list(self.__getitem__(roi,stranded=stranded))

    def __getitem__(self,roi,stranded=True):
        """Return list of features that overlap the region of interest (roi).
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            If `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
            Features that overlap `roi`

        Raises
        ------
        TypeError
            if `roi` is not a |GenomicSegment| or |SegmentChain|
        """
        if isinstance(roi,GenomicSegment):
            #roi_chain = SegmentChain(roi)
            roi_seg = roi
            roi_chain = SegmentChain(roi)
        elif isinstance(roi,SegmentChain):
            roi_chain = roi
            roi_seg = roi.spanning_segment
        else:
            raise TypeError("Query feature must be a GenomicSegment or SegmentChain")
        
        chrom = roi_seg.chrom
        feature_text = "\n".join(["\n".join(list(R.fetch(chrom,
                                                         X.start,
                                                         X.end))) \
                                                         for R in self.tabix_readers \
                                                         for X in roi_chain])
            
        features = (self._reader_class(cStringIO.StringIO(feature_text)))
        if stranded == True:
            features = [X for X in features if roi_chain.overlaps(X)]
        else:
            features = [X for X in features if roi_chain.unstranded_overlaps(X)]
        return features
