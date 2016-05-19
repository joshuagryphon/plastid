#!/usr/bin/env python
"""Constants, functions, and classes used by multiple readers in this subpackage


Functions & classes
-------------------
:func:`get_identical_attributes`
    Return a dictionary of all key-value pairs that are common to all `attr`
    dictionaries in a set of |SegmentChains|
    
|AssembledFeatureReader|
    Base class for readers that assemble high-level features (e.g. gapped
    alignments or transcripts) from one or more sub-features in an annotation file

"""
import itertools
import pysam
from plastid.util.io.filters import AbstractReader
from plastid.util.io.openers import NullWriter, multiopen
from plastid.genomics.roitools import GenomicSegment, SegmentChain, Transcript, add_three_for_stop_codon
from abc import abstractmethod


#===============================================================================
# INDEX: helper functions
#===============================================================================

def get_identical_attributes(features,exclude=None):
    """Return a dictionary of all key-value pairs that are identical for all |SegmentChains| in `features`
    
    Parameters
    ----------
    features : list 
        list of |SegmentChains|

    exclude : set
        attributes to exclude from identity criteria
    
    Returns
    -------
    dict
        Dictionary of all key-value pairs that have identical values in all the
        `attr` dictionaries of all the features in `features`
    """
    exclude = [] if exclude is None else exclude
    common_keys = set(features[0].attr.keys())
    for feature in features:
        common_keys &= set(feature.attr.keys())
    
    common_keys -= set(exclude)
    
    dtmp = { K : features[0].attr[K] for K in common_keys \
                 if all([X.attr[K] == features[0].attr[K] for X in features]) == True }
    return dtmp

#===============================================================================
# INDEX: classes
#===============================================================================

class AssembledFeatureReader(AbstractReader):
    """
    AssembledFeatureReader(*streams, return_type=SegmentChain, add_three_for_stop=False, tabix=False, printer=None, **kwargs)
    
    Abstract base class for readers that yield complex or discontinuous features
    such as transcripts or gapped alignments.
    
    For memory efficiency, all readers function as iterators. Readers built
    by subclassing |AssembledFeatureReader| are responsible for:
    
      - choosing when to yield assembled features
      
      - deciding how many subfeatures to hold in memory
      
      - overloading :meth:`~AssembledFeatureReader._assemble`
      

    Parameters
    ----------
    *streams : file-like
        One or more open filehandles of input data.
    
    return_type : |SegmentChain| or subclass, optional
        Type of feature to return from assembled subfeatures (Default: |SegmentChain|)

    add_three_for_stop : bool, optional
        Some annotation files exclude the stop codon from CDS annotations. If set to
        `True`, three nucleotides will be added to the threeprime end of each
        CDS annotation, **UNLESS** the annotated transcript contains explicit stop_codon 
        feature. (Default: `False`)
                    
    printer : file-like, optional
        Logger implementing a ``write()`` method. Default: |NullWriter|
    
    tabix : boolean, optional
        `streams` are `tabix`_-compressed (Default: `False`)

    **kwargs
        Other keyword arguments used by specific parsers
      
    
    Attributes
    ----------
    streams : file-like
        Input streams, usually constructed from or more open filehandles
    
    metadata : dict
        Various attributes gleaned from the stream, if any

    counter : int
        Cumulative line number counter over all streams

    printer : file-like, optional
        Logger implementing a ``write()`` method.

    return_type : class
        The type of object assembled by the reader. Typically an |SegmentChain|
        or a subclass thereof.
    
    rejected : list
        A list of transcript IDs that failed to assemble properly
    """
    
    def __init__(self,*streams,**kwargs):
        """
        AssembledFeatureReader(*streams, return_type=SegmentChain, add_three_for_stop=False, printer=None, tabix=False, **kwargs)
        
        Parameters
        ----------
        streams : file-like
            One or more fileneames or open filehandles of input data.
        
        return_type : |SegmentChain| or subclass, optional
            Type of feature to return from assembled subfeatures (Default: |SegmentChain|)

        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations. If set to
            `True`, three nucleotides will be added to the threeprime end of each
            CDS annotation, **UNLESS** the annotated transcript contains explicit stop_codon 
            feature. (Default: `False`)
                        
        printer : file-like, optional
            Logger implementing a ``write()`` method. Default: |NullWriter|
        
        tabix : boolean, optional
            `streams` point to `tabix`_-compressed files or are open
            :class:`~pysam.ctabix.tabix_file_iterator` (Default: `False`)

        **kwargs
            Other keyword arguments used by specific parsers
        """
        streams = multiopen(streams,fn=open,kwargs=dict(mode="rb"))
        
        if kwargs.get("tabix",False) == True:
            self.stream = itertools.chain.from_iterable((_tabix_iteradaptor(X) for X in streams))
        else:
            self.stream = itertools.chain.from_iterable(streams)

        self.counter = 0

        self.printer = kwargs.get("printer",NullWriter())

        self.return_type   = kwargs.get("return_type",SegmentChain)        
        add_three_for_stop = kwargs.get("add_three_for_stop",False)
        self._finalize =  add_three_for_stop_codon if add_three_for_stop == True else lambda x: x

        self.metadata    = {}
        self.rejected    = []
    
    @abstractmethod
    def _assemble(self,data):
        """Assemble features from data. This must be implemented in subclass.
        
        Returns
        -------
        |SegmentChain| or subclass
            Next feature assembled from `self.streams`, type specified by `self.return_type`
        """
    
    def filter(self,data):
        """Return next assembled feature from `self.stream`
        
        Returns
        -------
        |SegmentChain| or subclass
            Next feature assembled from `self.streams`, type specified by `self.return_type`
        """
        return self._finalize(self._assemble(data))


def _tabix_iteradaptor(stream):
    """Open `stream` as an iterator over a `tabix`_ file, returning raw strings from tabix data.
    
    Parameters
    ----------
    streams : open file-like, :class:`pysam.ctabix.tabix_file_iterator`
    
    Returns
    -------
    generator
        Generator of tab-delimited string records in `tabix`_ file
    """
    if not isinstance(stream,(pysam.ctabix.tabix_generic_iterator,
                              pysam.ctabix.tabix_file_iterator)
                      ):
        stream = pysam.tabix_file_iterator(stream,pysam.asTuple())
    
    return (str(X) for X in stream)
