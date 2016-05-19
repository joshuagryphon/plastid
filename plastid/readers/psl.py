#!/usr/bin/env python
"""This module defines a two classes for reading `PSL`_ files (made by, for example,
`blat`_):


|PSL_Reader|
    Read a `PSL`_ file line-by-line, converting each line into a |SegmentChain|
    or |Transcript|

|BundledPSL_Reader|
    Read `PSL`_ files, returning lists of |SegmentChains| grouped by query sequence.
"""
__date__ = "2011-09-01"
__author__ = "joshua"
import itertools
from plastid.readers.common import AssembledFeatureReader
from plastid.genomics.roitools import SegmentChain
from plastid.util.services.exceptions import FileFormatWarning, warn

class PSL_Reader(AssembledFeatureReader): 
    """
    PSL_reader(*streams, return_type=SegmentChain, add_three_for_stop=False, tabix=False, printer=None, **kwargs)
    
    Read `PSL`_ files into |SegmentChain|  or |Transcript| objects


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

    tabix : bool, optional
        `streams` point to `tabix`_-compressed files or are open
        :class:`~pysam.ctabix.tabix_file_iterator` (Default: `False`)


    **kwargs
        Other keyword arguments used by specific parsers

    
    Attributes
    ----------
    streams : file-like
        One or more open streams (usually filehandles) of input data.
    
    return_type : class
        The type of object assembled by the reader. Typically an |SegmentChain|
        or a subclass thereof. Must import a method called ``from_psl``

    counter : int
        Cumulative line number counter over all streams
    
    rejected : list
        A list of lines from `PSL`_ file that did not assemble properly
    
    metadata : dict
        Various attributes gleaned from the stream, if any    
    """
    def _assemble(self,line):
        """Read `PSL`_ files line-by-line into types specified by ``self.return_type``"""
        self.counter += 1
        if line.strip() == "":
            return self.__next()
        elif line.startswith("psLayout"):
            return self.__next__()
        elif line.lstrip().startswith("match"):
            return self.__next__()
        elif line.startswith("--"):
            return self.__next__()
        elif line.startswith("#"):
            return self.__next__()        
        else:
            try:
                return self.return_type.from_psl(line)
            except:
                self.rejected.append(line)
                warn("Rejecting line %s: %s" % (self.counter,line),FileFormatWarning)
                return self.__next__()        


class BundledPSL_Reader(PSL_Reader):
    """
    BundledPSL_reader(*streams, return_type=SegmentChain, add_three_for_stop=False, tabix=False, printer=None, **kwargs)
    
    Read `PSL`_ files, returning lists of |SegmentChains| grouped by query sequence.
    Use this when a given query sequence has multiple hits in your `PSL`_ file,
    and you want the output to be grouped.
    

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
    
    tabix : bool, optional
        `streams` point to `tabix`_-compressed files or are open
        :class:`~pysam.ctabix.tabix_file_iterator` (Default: `False`)

    **kwargs
        Other keyword arguments used by specific parsers    
    """

    def filter(self,line):
        """Process lines of `PSL`_ files input into |SegmentChain|, and group
        these by query sequence.
         
        Parameters
        ----------
        line : str
            line of `PSL`_ input
         
        Returns
        -------
        list
            list of |SegmentChain| objects sharing a query sequence 
        """
        ltmp = []
        aln = SegmentChain.from_psl(line)
        last_name = aln.attr["query_name"]
        try:
            while last_name == aln.attr["query_name"]:
                ltmp.append(aln)
                line = next(self.stream)
                aln = SegmentChain.from_psl(line)
                
            self.stream = itertools.chain([line],self.stream)
            return ltmp
        except StopIteration:
            # send final bundle
            return ltmp
