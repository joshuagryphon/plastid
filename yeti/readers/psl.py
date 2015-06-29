#!/usr/bin/env python
"""Contains tools for reading PSL alignment output from BLAT.
See http://pombe.nci.nih.gov/genome/goldenPath/help/blatSpec.html
for information on BLAT output 
"""
__date__ = "2011-09-01"
__author__ = "joshua"
from yeti.readers.common import AssembledFeatureReader
from yeti.genomics.roitools import SegmentChain
import itertools

class PSL_Reader(AssembledFeatureReader):
    """Reads BED files into |SegmentChain| or |Transcript|, saving metadata
    from track declarations et c. Malformed lines are stored in ``self.rejected``

    
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
        A list of transcript IDs that failed to assemble properly
    
    metadata : dict
        Various attributes gleaned from the stream, if any    
    """
    def _assemble(self,line):
        """Read PSL files line-by-line into types specified by ``self.return_type``"""
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
                self.printer.write("Rejecting line %s: %s" % (self.counter,line))
                return self.__next__()        


class BundledPSL_Reader(PSL_Reader):
    """Read PSL files, returning bundles of |SegmentChain| s grouped by query sequence.
    Use this when a given query sequence has multiple hits in your PSL file.
    """

    def filter(self,line):
        """Processes lines of BLAT input into |SegmentChain|, and groups
        these by query sequence.
         
        Parameters
        ----------
        line : str
            line of BLAT input to start parsing
         
        Returns
        -------
        list
            list of |SegmentChain| s for a particular query sequence
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