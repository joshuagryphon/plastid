#!/usr/bin/env python
"""This module contains |BED_Reader|, a parser that reads each line of a `BED`_
file into a |SegmentChain|, |Transcript|, or similar object. 

Examples
--------
Open a `BED`_ file, convert each line to a |Transcript|, and do something
with each transcript::

    >>> bed_reader = BED_Reader(open("some_file.bed"),return_type=Transcript)
    >>> for transcript in bed_reader:
            pass # do something fun

Retrieve a list of |SegmentChains| from a `BED`_ file::

    >>> my_chains = list(BED_Reader(open("some_file.bed"),return_type=SegmentChain))
    >>> my_chains[:5]
        [list of segment chains as output...]


Open several `Tabix`_-compressed `BED`_ files, and iterate over them as if
they were one stream::

    >>> import pysam
    >>> bed_files = [pysam.tabix_iterator(open(X), pysam.asTuple()) for X in ["file1.bed","file2.bed","file3.bed"]]
    >>> bed_reader = BED_Reader(*bed_files,tabix=True)
    >>> for segchain in bed_reader:
            pass # do something more interesting
                                

See Also
--------
`UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_.
    BED format specification at UCSC
"""
__date__ =  "Aug 23, 2011"
__author__ = "joshua"

import shlex
from yeti.readers.common import AssembledFeatureReader
from yeti.genomics.roitools import SegmentChain, Transcript
from yeti.util.services.decorators import deprecated, skipdoc

@skipdoc
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

@skipdoc
@deprecated
def BED_to_SegmentChain(stream):
    """Reads `BED`_ files line by line into |SegmentChain| objects
    
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
    """Reads `BED`_ files line-by-line into |SegmentChains| or |Transcripts|. 
    Metadata, if present in a track declaration, is saved in `self.metadata`.
    Malformed lines are stored in `self.rejected`, while parsing continues.

    
    Attributes
    ----------
    streams : file-like
        One or more open streams (usually filehandles) of input data.
    
    return_type : class
        The type of object assembled by the reader. Typically a |SegmentChain|
        or a subclass thereof. Must import a method called ``from_bed()``

    counter : int
        Cumulative line number counter over all streams
    
    rejected : list
        List of `BED`_ lines that could not be parsed
    
    metadata : dict
        Attributes declared in track line, if any
    """

    def _parse_track_line(self,inp):
        """Parse track line from `BED`_ file
        
        Parameters
        ----------
        inp : str
            track definition line from `BED`_ file
        
        Returns
        -------
        dict
            key-value pairs from `BED`_ line
        """
        ltmp = shlex.split(inp)
        for item in ltmp:
            k,v = item.split("=")
            self.metadata[k] = v
        
    def _assemble(self,line):
        """Read `BED`_ files line-by-line into types specified by `self.return_type`"""
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
