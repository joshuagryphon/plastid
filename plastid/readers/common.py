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
import copy
import itertools
from plastid.util.io.filters import AbstractReader
from plastid.util.io.openers import NullWriter
from plastid.genomics.roitools import GenomicSegment, SegmentChain, Transcript, add_three_for_stop_codon
from abc import abstractmethod


#===============================================================================
# INDEX: helper functions
#===============================================================================

def get_identical_attributes(features,exclude=set()):
    """Return a dictionary of all key-value pairs that are identical for all |SegmentChains| in `features`
    
    Parameters
    ----------
    features : list 
        list of |SegmentChains|
    
    Returns
    -------
    dict
        Dictionary of all key-value pairs that have identical values in all the
        `attr` dictionaries of all the features in `features`
    """
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
    """Abstract base class for readers that yield complex or discontinuous features
    such as transcripts or gapped alignments.
    
    For memory efficiency, all readers function as iterators. Readers built
    by subclassing |AssembledFeatureReader| are responsible for:
    
      - choosing when to yield assembled features
      
      - deciding how many subfeatures to hold in memory
      
      - overloading :meth:`~AssembledFeatureReader._assemble`
      
    
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
        Parameters
        ----------
        streams : file-like
            One or more open filehandles of input data.
        
        return_type : |SegmentChain| or subclass, optional
            Type of feature to return from assembled subfeatures (Default: |SegmentChain|)

        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations.
            If set to `True`, three nucleotides will be added to the threeprime
            end of each CDS annotation. (Default: `False`)
        
        printer : file-like, optional
            Logger implementing a ``write()`` method. Default: |NullWriter|
        
        tabix : boolean, optional
            `streams` are `tabix`_-compressed, and using the parser
            :py:class:`pysam.asTuple` (Default: `False`)

        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations. If set to
            `True`, three nucleotides will be added to the threeprime end of each
            CDS annotation, **UNLESS** the annotated transcript contains explicit stop_codon 
            feature. (Default: `False`)
                        
        **kwargs
            Other keyword arguments used by specific parsers
        """
        # this is a hack, because tabix iterator can no longer spit out
        # unparsed text
        self.stream = itertools.chain.from_iterable(streams)
        if kwargs.get("tabix",False) == True:
            self.stream = ("\t".join(X) for X in self.stream)

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
