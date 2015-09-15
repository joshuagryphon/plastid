"""Tools for looking up features in a region of interest within the a genome.

It is frequently useful to retrieve features that overlap specific regions 
of interest in the genome, for example, to find transcripts that overlap one
another. However, it would be inefficient to compare features that are
too far apart in the genome to overlap in the first place, or to scan through
a large annotation every time.

|GenomeHash|, |BigBedGenomeHash|, and |TabixGenomeHash| index annotations
by genomic location to avoid making unnecessary computations.

A |GenomeHash| may be created from a list or dictionary of features (e.g. |SegmentChains| or
|Transcripts|), or directly loaded from a genome annotation (in `BED`_, `GTF2`_,
`GFF3`_, or `PSL`_ format).

|BigBedGenomeHash| and |TabixGenomeHash| are more memory efficient
than |GenomeHash|, and take advantage of the indices already present in `BigBed`_ or 
`Tabix <http://samtools.sourceforge.net/tabix.shtml>`_-compressed files. They
thus avoid loading annotations into memory unless they are used.


Types of GenomeHash
-------------------
|GenomeHash|
    GenomeHash that can be created from memory-resident objects, such as a
    list or dict of |SegmentChains|, or an annotation file in `BED`_, `GTF2`_, `GFF3`_,
    or `PSL`_ format. Features are indexed by position in genome, allowing
    efficient lookup when, for example, comparing two genome annotations or
    comparing transcript isoforms.

|BigBedGenomeHash|
    A memory-efficient GenomeHash for use with `BigBed`_ files.

|TabixGenomeHash|
    A memory-efficient GenomeHash for use with  `BED`_, `GTF2`_, `GFF3`_, or `PSL`_
    files that have been compressed and indexed with `Tabix <http://samtools.sourceforge.net/tabix.shtml>`_.


Examples
--------
Create a |GenomeHash|::

    >>
    >>


Retrieve features overlapping a region of interest on chromosome II::

    >>
    >>
    
"""

#===============================================================================
# Memory-efficient ways to hash features across a genome
#===============================================================================
import copy
from plastid.util.services.mini2to3 import cStringIO
from plastid.util.io.openers import NullWriter
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
    
    Notes
    -----
    Because all features are stored in memory, for large genomes, a |TabixGenomeHash|
    or |BigBedGenomeHash| is much more memory-efficient. 
    """
    def __init__(self,features=[],binsize=DEFAULT_BIN_SIZE,do_copy=False):
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
    """A GenomeHash for `BigBed`_ files.
    
    Notes
    -----
    Because the underlying features are stored on disk, their data is immutable,
    and items cannot be added to the |BigBedGenomeHash| after it is initialized.
    
    Attributes
    ----------
    bigbedreader : |BigBedReader|
       |BigBedReader| connecting to BigBed file 
    """
    
    def __init__(self,filename,base_record_format="III",return_type=SegmentChain,cache_depth=5):
        """Create a |BigBedGenomeHash|
        
        Parameters
        ----------
        filename : filename
            Path to `BigBed`_ file (NOT open filehandle)

        base_record_format : str
            Format string for :py:func:`struct.unpack`, excluding endian-ness prefix
            and any notion of a null-terminated string. (Default: `III`)

        return_type : class implementing a :py:meth:`from_bed` method
            Class of object to return (Default: |SegmentChain|)
            
        cache_depth : int, optional
            Number of previously-fetched datablocks from `BigBed`_ file to keep
            resident in memory, to save time over repeated fetches to the same
            genomic regions. Increasing this number increases speed at the
            cost of increased memory use. (Default: `5`)
        """
        from plastid.readers.bigbed import BigBedReader
        self.filename = filename
        self.bigbedreader = BigBedReader(filename,
                                         base_record_format=base_record_format,
                                         return_type=return_type,
                                         cache_depth=cache_depth)
    
    def __del__(self):
        self.bigbedreader.close()
            
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
        return list(self.bigbedreader.__getitem__(roi,stranded=stranded))

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
        return list(self.bigbedreader[roi])
    

class TabixGenomeHash(AbstractGenomeHash):
    """A GenomeHash for a `tabix`_-indexed files
    
    Notes
    -----
    Because the underlying features are stored on disk, their data is immutable,
    and items cannot be added to the |TabixGenomeHash| after it is initialized.
    
    Attributes
    ----------
    filenames : list
        List of files used by |TabixGenomeHash|
        
    tabix_readers : list of :py:class:`pysam.Tabixfile`
       `Pysam`_ interfaces to underlying data files 
    """
    
    _READERS = { "GTF2" : GTF2_Reader,
                 "GFF3" : GFF3_Reader,
                 "BED"  : BED_Reader,
                 "PSL"  : PSL_Reader,
                }
    
    def __init__(self,filenames,data_format,printer=NullWriter()):
        """Create a |BigBedGenomeHash|
        
        Parameters
        ----------
        filenames : list of str
            One or more paths to `Tabix`_-compressed files (NOT open filehandle)

        data : str
            Format of tabix-compressed file(s). Choices are:
            `'GTF2'`,`'GFF3'`,`'BED'`,`'PSL'`
        """
        from pysam import Tabixfile
        self.filenames = filenames
        self.printer = printer
        try:
            self._reader_class = TabixGenomeHash._READERS[data_format]
        except ValueError:
            msg = "Supported file formats for TabixGenomeHash are: %s" % ", ".join(sorted(TabixGenomeHash._READERS.keys()))
            printer.write(msg)
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
