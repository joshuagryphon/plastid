#!/usr/bin/env python
"""This module contains |BED_Reader|, an iterator that reads each line of a `BED`_
or :term:`extended BED` file into a |SegmentChain|, |Transcript|, or similar object. 

.. contents::
   :local:

Module contents
---------------

.. autosummary::

   BED_Reader
   bed_x_formats

Examples
--------
Read entries in a `BED`_ file as |Transcripts|. `thickEnd` and `thickStart`
columns will be interpreted as the endpoints of coding regions::

    >>> bed_reader = BED_Reader("some_file.bed",return_type=Transcript)
    >>> for transcript in bed_reader:
            pass # do something fun with each Transcript/SegmentChain

If `return_type` is unspecified, `BED`_ lines are read as |SegmentChains|::

    >>> my_chains = list(BED_Reader("some_file.bed"))
    >>> my_chains[:5]
        [list of segment chains as output...]

Open an :term:`extended BED` file, which contains additional columns for `gene_id`
and `favorite_color`. Values for these attributes will be stored in the `attr`
dict of each |Transcript|::

    >>> bed_reader = BED_Reader("some_file.bed",return_type=Transcript,extra_columns=["gene_id","favorite_color"])

Open several `Tabix`_-compressed `BED`_ files, and iterate over them as if
they were one stream::

    >>> bed_reader = BED_Reader("file1.bed.gz","file2.bed.gz",tabix=True)
    >>> for chain in bed_reader:
    >>>     pass # do something interesting with each chain
                                

See Also
--------
`UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_.
    BED format specification at UCSC
"""
__date__ =  "Aug 23, 2011"
__author__ = "joshua"

import shlex
from plastid.readers.common import AssembledFeatureReader
from plastid.util.services.exceptions import FileFormatWarning, warn


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
"""Column names and types for various :term:`extended BED` formats used by 
the `ENCODE`_ project. These can be passed to the `extra_columns` keyword of
:class:`BED_Reader`.""" 


class BED_Reader(AssembledFeatureReader):
    """
    BED_Reader(*streams, return_type=SegmentChain, add_three_for_stop=False, extra_columns=0, printer=None, tabix=False)
    
    Reads `BED`_ and :term:`extended BED` files line-by-line into |SegmentChains|
    or |Transcripts|. Metadata, if present in a track declaration, is saved
    in `self.metadata`. Malformed lines are stored in `self.rejected`, while
    parsing continues.

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

    extra_columns: int or list optional
        Extra, non-BED columns in :term:`extended BED` format file corresponding
        to feature attributes. This is common in `ENCODE`_-specific `BED`_ variants.
        
        if `extra-columns` is:
        
          - an :class:`int`: it is taken to be the
            number of attribute columns. Attributes will be stored in
            the `attr` dictionary of the |SegmentChain|, under names like
            `custom0`, `custom1`, ... , `customN`.

          - a :class:`list` of :class:`str`, it is taken to be the names
            of the attribute columns, in order, from left to right in the file.
            In this case, attributes in extra columns will be stored under
            their respective names in the `attr` dict.

          - a :class:`list` of :class:`tuple`, each tuple is taken
            to be a pair of `(attribute_name, formatter_func)`. In this case,
            the value of `attribute_name` in the `attr` dict of the |SegmentChain|
            will be set to `formatter_func(column_value)`.
        
        (Default: 0)
                    
    printer : file-like, optional
        Logger implementing a ``write()`` method. Default: |NullWriter|
    
    tabix : boolean, optional
        `streams` point to `tabix`_-compressed files or are open
        :class:`~pysam.ctabix.tabix_file_iterator` (Default: `False`)
        

    Examples
    --------
    Read entries in a `BED`_ file as |Transcripts|. `thickEnd` and `thickStart`
    columns will be interpreted as the endpoints of coding regions::
    
        >>> bed_reader = BED_Reader(open("some_file.bed"),return_type=Transcript)
        >>> for transcript in bed_reader:
        >>>     pass # do something fun
    
    Open an :term:`extended BED` file that contains additional columns for `gene_id`
    and `favorite_color`. Values for these attributes will be stored in the `attr`
    dict of each |Transcript|::
    
        >>> bed_reader = BED_Reader(open("some_file.bed"),return_type=Transcript,extra_columns=["gene_id","favorite_color"])

    Open several `Tabix`_-compressed `BED`_ files, and iterate over them as if
    they were one uncompressed stream::
    
        >>> bed_reader = BED_Reader("file1.bed.gz","file2.bed.gz",tabix=True)
        >>> for chain in bed_reader:
        >>>     pass # do something more interesting

    
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
        Extra, non-BED columns in :term:`extended BED` format file corresponding to feature
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
            types presently include:

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
        """
        BED_Reader(*streams, return_type=SegmentChain, add_three_for_stop=False, extra_columns=0, printer=None, tabix=False)
        
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

        extra_columns: int or list optional
            Extra, non-BED columns in :term:`Extended BED`_ format file corresponding
            to feature attributes. This is common in `ENCODE`_-specific `BED`_ variants.
            
            if `extra-columns` is:
            
              - an :class:`int`: it is taken to be the
                number of attribute columns. Attributes will be stored in
                the `attr` dictionary of the |SegmentChain|, under names like
                `custom0`, `custom1`, ... , `customN`.

              - a :class:`list` of :class:`str`, it is taken to be the names
                of the attribute columns, in order, from left to right in the file.
                In this case, attributes in extra columns will be stored under
                their respective names in the `attr` dict.

              - a :class:`list` of :class:`tuple`, each tuple is taken
                to be a pair of `(attribute_name, formatter_func)`. In this case,
                the value of `attribute_name` in the `attr` dict of the |SegmentChain|
                will be set to `formatter_func(column_value)`.
            
            (Default: 0)
                        
        printer : file-like, optional
            Logger implementing a ``write()`` method. Default: |NullWriter|
        
        tabix : boolean, optional
            `streams` are `tabix`_-compressed (Default: `False`)
        """
        AssembledFeatureReader.__init__(self,*args,**kwargs)
        self.extra_columns = kwargs.get("extra_columns",0)

    def _parse_track_line(self,inp):
        """Parse track line from `BED`_ / extended BED file
        
        Parameters
        ----------
        inp : str
            track definition line from `BED`_  / extended BED file

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
                    warn("Extra columns specified by %s track type declaration (%s) don't match those specified by user (%s). Using those specified by user." %\
                         (track_type,track_format_columns,my_columns),FileFormatWarning)
                    self.metadata["type"] = "custom"
            else:
                self.printer.write("Found track type '%s' in track definition line." % track_type)
        
    def _get_extra_column_names(self):
        """Return names of extra columns in extended BED file)"""
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
                    msg += "Maybe this BED has extra columns (i.e. is an extended BED file)?"

                msg += ("\n    %s" % line)
                warn(msg,FileFormatWarning)
                return self.__next__()
