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
import warnings
from plastid.readers.common import AssembledFeatureReader
from plastid.genomics.roitools import SegmentChain, Transcript
from plastid.util.services.decorators import deprecated, skipdoc
from plastid.util.services.exceptions import FileFormatWarning

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

bed_x_formats = {
    "bedDetail" : [("ID",str),
                   ("description",str),
                  ],
    "narrowPeak" : [("signalValue",float),
                    ("pValue",float),
                    ("qValue",float),
                    ("peak",int)],
    "broadPeak"  : [("signalValue",float),
                    ("pValue",float),
                    ("qValue",float)],
    "gappedPeak" : [("signalValue",float),
                    ("pValue",float),
                    ("qValue",float)],
    "tagAlign"   : [("sequence",str),
                    ("score",float),
                    ("strand",str)],
    "pairedTagAlign" : [("seq1",str),
                        ("seq2",str)],
    "peptideMapping" : [("rawScore",float),
                        ("spectrumId",str),
                        ("peptideRank",int),
                        ("peptideRepeatCount",int)],

}

class BED_Reader(AssembledFeatureReader):
    """Reads `BED`_ and :term:`BED X+Y` files line-by-line into |SegmentChains|
    or |Transcripts|. Metadata, if present in a track declaration, is saved
    in `self.metadata`. Malformed lines are stored in `self.rejected`, while
    parsing continues.

    
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
    
    extra_columns : int or list, optional
        Extra, non-BED columns in :term:`BED X+Y`_ format file corresponding to feature
        attributes. This is common in `ENCODE`_-specific `BED`_ variants.
        
        if `extra_columns` is:
        
          - an :class:`int`: it is taken to be the
            number of attribute columns. Attributes will be stored in
            the `attr` dictionary of the |SegmentChain|, under names like
            `custom0`, `custom1`, ... , `customN`.

          - a :class:`list` of :class:`str`, it is taken to be the names
            of the attribute columns, in order, from left to right in the file.
            In this case, attributes in extra columns will be stored under
            there respective names in the `attr` dict.

          - a :class:`list` of :class:`tuple`, each tuple is taken
            to be a pair of `(attribute_name, formatter_func)`. In this case,
            the value of `attribute_name` in the `attr` dict of the |SegmentChain|
            will be set to `formatter_func(column_value)`.
        
        If unspecified, :class:`BED_Reader` reads the track declaration line
        (if present), and:

          - if a known track type is specified by the `type` field, it attempts
            to format the extra columns as specified by that type. Known track
            types prently include:

              - bedDetail
              - narrowPeak
              - broadPeak
              - gappedPeak
              - tagAlign
              - pairedTagAlign
              - peptideMapping
          
          - if not, it assumes 0 non-`BED`_ fields are present, and that all columns
            are `BED`_ formatted.
    """
    def __init__(self,*args,**kwargs):
        AssembledFeatureReader.__init__(self,*args,**kwargs)
        self.extra_columns = kwargs.get("extra_columns",0)

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
        self.metadata = {}
        ltmp = shlex.split(inp.strip("\n"))
        for item in ltmp:
            k,v = item.split("=")
            self.metadata[k] = v

        track_type = self.metadata.get("type",None)
        if track_type is not None:
            if track_type in bed_x_formats:
                self.printer.write("Found track type '%s' in track definition line. Assuming extra columns follow UCSC definitions." % track_type)
                if self.extra_columns == 0:
                    self.extra_columns = bed_x_formats[track_type]
                elif self.extra_columns != bed_x_formats[track_type]:
                    my_columns = self._get_extra_column_names()
                    track_format_columns = ",".join([X[0] for X in bed_x_formats[track_type]])
                    warnings.warn("Extra columns specified by %s track type declaration (%s) don't match those specified by user (%s). Using those specified by user." %\
                                  (track_type,track_format_columns,my_columns),FileFormatWarning)
                    self.metadata["type"] = "custom"
            else:
                self.printer.write("Found track type '%s' in track definition line." % track_type)
        
    def _get_extra_column_names(self):
        """Return names of extra columns in BED X+Y file)"""
        if isinstance(self.extra_columns,int):
            my_columns = "%s unnamed columns" % self.extra_columns
        elif isinstance(self.extra_columns,list):
            if all([isinstance(X,tuple) for X in self.extra_columns]):
                my_columns = ",".join([X[0] for X in self.extra_columns])
            elif all([isinstance(X,str) for X in self.extra_columns]):
                my_columns = ",".join(self.extra_columns)

        return my_columns

    def _assemble(self,line):
        """Read `BED`_ files line-by-line into types specified by `self.return_type`"""
        self.counter += 1
        if line.strip() == "":
            return self.__next__()
        elif line.startswith("browser"):
            return self.__next__()
        elif line.startswith("track"):
            # reset metadata
            self._parse_track_line(line[5:])
            return self.__next__()
        elif line.startswith("#"):
            return self.__next__()
        else:
            try:
                return self.return_type.from_bed(line,extra_columns=self.extra_columns)
            except:
                self.rejected.append(line)
                msg = "Cannot parse BED line number %s. " % self.counter
                if self.metadata.get("type",None) is not None:
                    msg += ("Are you sure this is a %s BED file with extra columns (%s)?" % (self.metadata.get("type"),self._get_extra_column_names()))
                elif self.extra_columns != 0:
                    msg += ("Are you sure this BED file has extra columns (%s)?" % self._get_extra_column_names())
                else:
                    msg += "Maybe this BED has extra columns (i.e. is a BED X+Y file)?"

                msg += ("\n    %s" % line)
                warnings.warn(msg,FileFormatWarning)
                return self.__next__()
