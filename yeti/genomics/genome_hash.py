"""Tools for quickly locating features of interest in a sub-region of a genome.

It is frequently useful to retrieve features that overlap specific regions 
of interest in the genome, for example, to find transcripts that overlap one
another. However, it would be inefficient to compare features that are
too far apart in the genome to overlap in the first place. 

|GenomeHash|, |BigBedGenomeHash|, and |TabixGenomeHash| index annotations
by genomic location to avoid making unnecessary comparisons. A |GenomeHash|
may be created from a list or dictionary of features (e.g. |SegmentChain| s or
|Transcript| s), or directly loaded from a genome annotation (in `BED`_, `GTF2`_,
`GFF3`_, or `PSL`_ format).

|BigBedGenomeHash| and |TabixGenomeHash| objects are more memory efficient,
and take advantage of the indices already present in `BigBed`_ or Tabix
<http://samtools.sourceforge.net/tabix.shtml>`_-compressed files, to avoid
loading annotations into memory before they are used (if they even are at all).


Important classes
-----------------
|GenomeHash|
    GenomeHash that can be created from memory-resident objects, such as a
    list or dict of |SegmentChain| s, or an annotation file in `BED`_, `GTF2`_, `GFF3`_,
    or `PSL`_ format. Indexes regions of interest by position in genome, allowing
    efficient lookup of nearby or overlapping features when, for example,
    comparing two genome annotations or comparing transcript isoforms

|BigBedGenomeHash|
    A more memory-efficient implementation of |GenomeHash| for use with `BigBed`_ files,
    which contain their own genomic indices. Therefore, features only need to be
    loaded into memory if and when they are be used.

|TabixGenomeHash|
    A more memory-efficient implementation of |GenomeHash| for use with 
    `BED`_, `GTF2`_, `GFF3`_, or `PSL`_ files that have been compressed and
    indexed with `Tabix <http://samtools.sourceforge.net/tabix.shtml>`_.


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
#import cStringIO
from yeti.util.services.mini2to3 import cStringIO
from pysam import Tabixfile
from yeti.util.io.openers import NullWriter
from yeti.readers.bed import BED_Reader
from yeti.readers.gff import GTF2_Reader, GFF3_Reader
from yeti.readers.psl import PSL_Reader
from yeti.genomics.roitools import GenomicSegment, SegmentChain
from yeti.readers.bigbed import BigBedReader
from abc import abstractmethod

DEFAULT_BIN_SIZE=20000


class AbstractGenomeHash(object):
    """Abstract base class for hashes that map |SegmentChain| objects to
    neighborhoods, allowing quick lookup for comparisons of overlap
    or other behavior
    """
    @abstractmethod
    def get_overlapping_features(self,feature,stranded=True):
        """Return list of features in all the bins occupied by 'feature'
        
        Parameters
        ----------
        feature : |GenomicSegment| or |SegmentChain|
            Query feature

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
            Overlapping |SegmentChain| s.

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        pass

    @abstractmethod
    def __getitem__(self,roi):
        """Return list of features that overlap the region of interest (roi),
        minding the strand.
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
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
    """Creates a hash mapping a list of dictionary of |SegmentChain| 
    to neighborhoods, allowing quick lookup for comparisons of overlap
    or other behavior.
    
    Notes
    -----
    All features are stored in memory. For large genomes, a |TabixGenomeGash|
    or |BigBedGenomeHash| is much more memory-efficient. 
    """
    def __init__(self,features,binsize=DEFAULT_BIN_SIZE,do_copy=True):
        """Create a |GenomeHash|
        
         Parameters
         ----------
         features : dict or list
             dict or list of |SegmentChain|
            
         binsize : int
             size of "neighborhood" for hash.
            
         do_copy : bool
             If ``True`` (default), the internal `feature_dict`
             used to make the genome_hash will be a copy
             of the one supplied in `feature_dict` parameter.                         
             This comes at a speed cost, but will prevent side effects
             if the features are changed outside the hash. 
             
             If ``False``, creation of the |GenomeHash| will be much faster.
        """
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
            dict or list of |SegmentChain|
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
        """Create the |GenomeHash|. Internally, the hash is of type
        dict[chrom][strand] = list<feature_name>
           
        Returns
        -------
        dict
            Hierarchichal dictionary: dict[chrom][strand] = list<str>
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
        
    def _get_hash_bins(self,feature):
        """Returns a list of genome bins that a given feature falls into
        
        Parameters
        ----------
        feature : |GenomicSegment| or |SegmentChain|
            Query feature
           
           
        Returns
        -------
        list<int>
            List of hash bins that feature covers
               
           
        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        if isinstance(feature,GenomicSegment):
            iv = feature
        elif isinstance(feature,SegmentChain):
            iv = feature.spanning_segment
        else:
            raise TypeError("Query feature must be a GenomicSegment or SegmentChain")
        bin_start = iv.start // self.binsize
        bin_end = iv.end // self.binsize
        return list(range(bin_start,bin_end+1))
    
    def _get_nearby_feature_ids(self,feature,stranded=True):
        """Return unique IDs of |SegmentChain| s in all the bins occupied by 'feature'
        
        Parameters
        ----------
        feature : |GenomicSegment| or |SegmentChain|
           Query feature
           
        stranded : bool
            If ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
           
                             
        Returns
        -------
        list<str>
               Unique IDs of features near (within ``self.binsize`` distance from) ``feature``
           
           
        Raises
        ------
        TypeError
            if feature is not |GenomicSegment| or |SegmentChain|
        """
        if isinstance(feature,GenomicSegment):
            iv = feature
        elif isinstance(feature,SegmentChain):
            iv = feature.spanning_segment
        else:
            raise TypeError("Query feature must be a GenomicSegment or SegmentChain")
        feature_ranges = self._get_hash_bins(feature)
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

    def get_nearby_features(self,feature,stranded=True):
        """Return list of features in all the bins occupied by 'feature'
        
        Parameters
        ----------
        feature : |GenomicSegment| or |SegmentChain|
            Query feature

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list<str>
               Names of features near (within ``self.binsize`` distance from) ``feature``
               

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        nearby_feature_ids = self._get_nearby_feature_ids(feature,stranded=stranded)
        return [self.feature_dict[X] for X in nearby_feature_ids]

    def get_nearby_feature_names(self,feature,stranded=True):
        """Return list of the names of features in all the bins occupied by 'feature'
        
        Parameters
        ----------
        feature : |GenomicSegment| or |SegmentChain|
            Query feature

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
            Features near (within ``self.binsize`` distance from) ``feature``

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        nearby_feature_ids = self._get_nearby_feature_ids(feature,stranded=stranded)
        return [self._id_to_names[X] for X in nearby_feature_ids]

    def get_overlapping_features(self,roi,stranded=True):
        """Return list of features in all the bins occupied by ``roi``
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
               Overlapping features

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        nearby_features = self.get_nearby_features(roi, stranded=stranded)
        if stranded == False:
            fn = roi.unstranded_overlaps
        else:
            fn = roi.overlaps
        return [X for X in nearby_features if fn(X) == True]
    
    def __getitem__(self,roi):
        """Return list of features that overlap the region of interest (roi),
        minding the strand.
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
               Overlapping features

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| |SegmentChain|
        """
        return self.get_overlapping_features(roi,stranded=True)
                

class BigBedGenomeHash(AbstractGenomeHash):
    """A fast, memory-efficient |AbstractGenomeHash| built upon a randomly-accessible BigBed file.
    
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
        
        filename : filename
            Path to BigBed file (NOT open filehandle)

        base_record_format : str
            Format string for :py:func:'struct.unpack`, excluding endian-ness prefix
            and any notion of a null-terminated string. (Default: "III")

        return_type : class implementing a :py:meth:`from_bed` method
            Class of object to return (Default: |SegmentChain|)
            
        cache_depth : int, optional
            Number of previously-fetched datablocks from BigBed file to keep
            resident in memory, to save time over repeated fetches to the same
            genomic regions. Increasing this number increases speed at the
            cost of memory. (Default: 5)
        """
        self.filename = filename
        self.bigbedreader = BigBedReader(filename,
                                         base_record_format=base_record_format,
                                         return_type=return_type,
                                         cache_depth=cache_depth)
    
    def __del__(self):
        self.bigbedreader.close()
            
    def get_overlapping_features(self,roi,stranded=True):
        """Return list of features overlapping ``roi``
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
            Overlapping features

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        return list(self.bigbedreader.__getitem__(roi,stranded=stranded))

    def __getitem__(self,roi,stranded=True):
        """Return list of features that overlap the region of interest (roi),
        minding the strand.
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
            Overlapping features

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        return list(self.bigbedreader[roi])
    

class TabixGenomeHash(AbstractGenomeHash):
    """A memory-efficient |AbstractGenomeHash| built upon a tabix-indexed file.
    
    Notes
    -----
    Because the underlying features are stored on disk, their data is immutable,
    and items cannot be added to the |TabixGenomeHash| after it is initialized.
    
    Attributes
    ----------
    filenames : list
        List of files used by |TabixGenomeHash|
        
    tabix_reader : :py:class:`Tabixfile`
       Pysam interface to underlying data 
    """
    
    _READERS = { "GTF2" : GTF2_Reader,
                 "GFF3" : GFF3_Reader,
                 "BED"  : BED_Reader,
                 "PSL"  : PSL_Reader,
                }
    
    def __init__(self,filenames,data_format,printer=NullWriter()):
        """Create a |BigBedGenomeHash|
        
        filenames : list of str
            One or more paths to Tabix files (NOT open filehandle)

        data : str
            Format of tabix-compressed file(s). Choices are:
            *GTF2*,*GFF3*,*BED*,*PSL*
        """
        self.filenames = filenames
        self.printer = printer
        try:
            self._reader_class = TabixGenomeHash._READERS[data_format]
        except KeyError:
            msg = "Supported formats for TabixGenomeHash are: %s" % ", ".join(sorted(TabixGenomeHash._READERS.keys()))
            printer.write(msg)
            raise KeyError(msg)
        
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
        """Return list of features overlapping ``roi``
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
        Returns
        -------
        list
               Overlapping features

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        return list(self.__getitem__(roi,stranded=stranded))

    def __getitem__(self,roi,stranded=True):
        """Return list of features that overlap the region of interest (roi),
        minding the strand.
        
        Parameters
        ----------
        roi : |GenomicSegment| or |SegmentChain|
            Query feature indicating region of interest

        stranded : bool
            if ``True``, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands
                             
                             
        Returns
        -------
        list
            Overlapping features

        Raises
        ------
        TypeError
            if feature is not a |GenomicSegment| or |SegmentChain|
        """
        if isinstance(roi,GenomicSegment):
            roi_ivc = SegmentChain(roi)
        elif isinstance(roi,SegmentChain):
            roi_ivc = roi
        else:
            raise TypeError("Query feature must be a GenomicSegment or SegmentChain")
        
        feature_text = "\n".join(["\n".join(list(R.fetch(roi_ivc.spanning_segment.chrom,
                                                         roi_ivc.spanning_segment.start,
                                                         roi_ivc.spanning_segment.end))) \
                                                         for R in self.tabix_readers])
            
        features = list(self._reader_class(cStringIO.StringIO(feature_text)))
        if stranded == True:
            features = [X for X in features if roi_ivc.overlaps(X)]
        else:
            features = [X for X in features if roi_ivc.unstranded_overlaps(X)]
        return features
