#!/usr/bin/env python
"""This module defines a two classes for reading `PSL`_ files (made by, for example,
`blat`_):


|PSL_Reader|
    Read a `PSL`_ file line-by-line, converting each line into a |SegmentChain|
    or |Transcript|

|BundledPSL_Reader|
    Read `PSL`_ files, returning lists of |SegmentChains| grouped by query sequence.


Examples
--------

Read individual entries from a `PSL`_ file into |SegmentChain| objects::

    >>>


Group multiple entries from a `PSL`_ file by query sequence::

    >>>


"""
__date__ = "2011-09-01"
__author__ = "joshua"
import warnings
import itertools
from plastid.readers.common import AssembledFeatureReader
from plastid.genomics.roitools import SegmentChain
from plastid.util.services.exceptions import FileFormatWarning

class PSL_Reader(AssembledFeatureReader): 
    """Read `PSL`_ files into |SegmentChain|  or |Transcript| objects
    
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
                warnings.warn("Rejecting line %s: %s" % (self.counter,line),FileFormatWarning)
                return self.__next__()        


class BundledPSL_Reader(PSL_Reader):
    """Read `PSL`_ files, returning lists of |SegmentChains| grouped by query sequence.
    Use this when a given query sequence has multiple hits in your `PSL`_ file,
    and you want the output to be grouped.
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
