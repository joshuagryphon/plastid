#!/usr/bin/env python
"""Constants, functions, and classes used by readers of multiple file types

Important methods & classes
---------------------------
|AssembledFeatureReader|
    Base class for readers that assemble high-level features (e.g. gapped
    alignments or transcripts) from one or more sub-features in an annotation file

"""
import copy
import itertools
from yeti.util.io.filters import AbstractReader
from yeti.util.io.openers import NullWriter
from yeti.genomics.roitools import SegmentChain, Transcript
from abc import abstractmethod


#===============================================================================
# INDEX: helper functions
#===============================================================================

def add_three_for_stop_codon(tx):
    """Extend an annotated CDS region, if present, by three nucleotides
    at the threeprime end. Use in cases when annotation files exclude the stop
    codon from the annotated CDS.
    
    Parameters
    ----------
    tx : |Transcript|
        query transcript
        
    Returns
    -------
    |Transcript|
        |Transcript| with same attributes as `tx`, but with a longer CDS
    
    Raises
    ------
    IndexError
        if a three prime UTR is defined that terminates before the complete stop codon
    """
    if tx.cds_genome_end is not None and tx.cds_genome_start is not None:
        segments  = copy.deepcopy(tx._segments)
        attr = copy.deepcopy(tx.attr)
        if tx.spanning_segment.strand == "+":
            # if cds at end of transcript, both need to grow by three
            if tx.spanning_segment.end == tx.cds_genome_end:
                segments[-1].end += 3
                attr["cds_genome_end"] = tx.cds_genome_end + 3
            else:
                # if new end will be end of transcript, set manually
                # because key won't be in position hash
                if tx.cds_end + 3 == tx.get_length():
                    new_end = tx.spanning_segment.end
                else:
                    # position is internal, so we fetch
                    new_end = tx.get_genomic_coordinate(tx.cds_end+3)[1]
                    # correct for rare cases in which stop codons and at the end
                    # of a given exon interval, and are falsely mapped to first position
                    # of next exon. Instead, we should map to 1 + the genomic position
                    # of the end of the given exon.
                    for n,segment in enumerate(tx):
                        if new_end == segment.start:
                            new_end = tx[n-1].end
                    
                attr["cds_genome_end"] = new_end
            return Transcript(*segments,**attr)
        else:
            # if cds starts at beginning of transcript, both need to grow by three
            if tx.spanning_segment.start == tx.cds_genome_start:
                segments[0].start -= 3
                attr["cds_genome_start"] = tx.cds_genome_start - 3
            else:
                # otherwise, fetch
                # minus one here might look funny, but is due to minus-strand
                # feature-ness
                new_start = tx.get_genomic_coordinate(tx.cds_end + 3 - 1)[1]
                attr["cds_genome_start"] = new_start
            return Transcript(*segments,**attr)
    else:
        return tx


def get_identical_attributes(features,exclude=set()):
    """Return a dictionary of all key-value pairs that are identical for all features in `features`
    
    Parameters
    ----------
    features : list 
        list of |SegmentChain| objects representing features to consider
    
    Returns
    -------
    dict
        Dictionary of all key-value pairs that have identical values in all
        features in `features`
    """
    common_keys = set(features[0].attr.keys())
    for feature in features:
        common_keys &= set(feature.attr.keys())
    
    common_keys -= set(exclude)
    
    dtmp = { K : features[0].attr[K] for K in common_keys \
                 if len(set([X.attr[K] for X in features])) == 1  }
    return dtmp

#===============================================================================
# INDEX: classes
#===============================================================================

class AssembledFeatureReader(AbstractReader):
    """Abstract base class for readers that assemble features in various
    annotation file formats into complex or discontinuous features, such
    as transcripts or gapped alignments. For memory efficiency, all readers
    function as iterators over their respective streams. Choosing when to
    process sub-features, and how many to hold in memory, is up to implementation
    in subclasses, which should additionally define `self._assemble`.
    
    Attributes
    ----------
    stream : file-like
        Input stream, usually constructed from or more open filehandles
    
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
            Set to *True* if incoming streams are `tabix`_-compressed, and
            using the parser :py:class:`pysam.asTuple` (Default: `False`)

        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations. If set to
            *True*, three nucleotides will be added to the threeprime end of each
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
            Next feature assembled from `self.stream`, specified by `self.return_type`
        """
    
    def filter(self,data):
        """Return next assembled feature from `self.stream`
        
        Returns
        -------
        |SegmentChain| or subclass
            Next feature assembled from `self.stream`, specified by `self.return_type`
        """
        return self._finalize(self._assemble(data))
