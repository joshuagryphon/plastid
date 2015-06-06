#!/usr/bin/env python
"""Utilities for reading BED files 

    :py:func:`BED_to_Transcripts`
        Read BED files to |Transcript| objects
    
    :py:func:`BED_to_SegmentChain`
        Reads BED files to |SegmentChain| objects
    
    |BED_Reader|
        Reads BED files to |SegmentChain| objects

See Also
--------
`UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_.
    BED specification at UCSC
"""
__date__ =  "Aug 23, 2011"
__author__ = "joshua"

import shlex
from yeti.readers.common import AssembledFeatureReader
from yeti.genomics.roitools import SegmentChain, Transcript
from yeti.util.services.decorators import deprecated

@deprecated
def BED_to_Transcripts(stream,add_three_for_stop=False):
    """Reads BED files line by line into |Transcript| objects
    
    Parameters
    ----------
    stream : file-like
        Stream of BED4-BED12 format data

    add_three_for_stop : bool, optional
        Some annotation files exclude the stop codon from CDS annotations. If set to
        True, three nucleotides will be added to the threeprime end of each
        CDS annotation. Default: False
    
    Yields
    ------
    |Transcript|
    """
    reader = BED_Reader(stream,return_type=Transcript,add_three_for_stop=add_three_for_stop)
    for ivc in reader:
        yield ivc

@deprecated
def BED_to_SegmentChain(stream):
    """Reads BED files line by line into |SegmentChain| objects
    
    Parameters
    ----------
    stream : file-like
        Stream of BED4-BED12 format data
    
    Yields
    ------
    |SegmentChain|
    """
    reader = BED_Reader(stream,return_type=SegmentChain)
    for ivc in reader:
        yield ivc

class BED_Reader(AssembledFeatureReader):
    """Reads BED files into |SegmentChain| or |Transcript|, saving metadata
    from track declarations et c. Malformed lines are stored in ``self.rejected``

    
    Attributes
    ----------
    streams : file-like
        One or more open streams (usually filehandles) of input data.
    
    return_type : class
        The type of object assembled by the reader. Typically an |SegmentChain|
        or a subclass thereof. Must import a method called ``from_bed``

    counter : int
        Cumulative line number counter over all streams
    
    rejected : list
        A list of transcript IDs that failed to assemble properly
    
    metadata : dict
        Attributes declared in track line, if any    
    """

    def _parse_track_line(self,inp):
        """Parse tokens in track declaration
        
        Parameters
        ----------
        inp : str
            track definition line from BED file
        
        Returns
        -------
        dict : key-value pairs from BED line
        """
        ltmp = shlex.split(inp)
        for item in ltmp:
            k,v = item.split("=")
            self.metadata[k] = v
        
    def _assemble(self,line):
        """Read BED files line-by-line into types specified by ``self.return_type``"""
        self.counter += 1
        if line.strip() == "":
            return self.__next__()
        elif line.startswith("browser"):
            return self.__next__()
        elif line.startswith("track"):
            self._parse_track_line(line[5:])
            return self.__next__()
        elif line.startswith("#"):
            return self.__next__()
        else:
            try:
                return self.return_type.from_bed(line)
            except:
                self.rejected.append(line)
                self.printer.write("Rejecting line %s: %s" % (self.counter,line))
                return self.__next__()
