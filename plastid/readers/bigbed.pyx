#!/usr/bin/env python
"""|BigBedReader|, a parser for `BigBed`_ files.

.. contents::
   :local:


Summary
-------

In contrast to `BED`_, `GTF2`_, and `GFF3`_ files, `BigBed`_ files are binary,
indexed, and randomly-accessible. This means:

    - |BigBedReader| can be used to iterate over records, like a reader, **or**
      to fetch records that cover a region of interest, in the manner of a
      |GenomeHash|

    - `BigBed`_ use less memory, because their records don't need to be loaded
      into memory to be parsed or accessed.

    - Indexes `BigBed`_ files can be searched for matching records


Module Contents
---------------

.. autosummary::

   BigBedReader
   BigBedIterator


Examples
--------
Iterate over all features in a `BigBed`_ file::

    >>> my_reader = BigBedReader("some_file.bb", return_type=Transcript)
    >>> for feature in my_reader:
    >>>    pass # do something with each Transcript

`BigBed`_ files can be accessed as dictionaries. To find features overlapping a
region of interest::

    >>> roi = GenomicSegment("chrI", 0, 100000, "+")
    >>> overlapping_features = my_reader[roi]
    >>> list(overlapping_features)
    [ list of SegmentChains/Transcripts ]

Find features that match keyword(s) in a certain field::

    >>> # which fields are indexed and searchable?
    >>> my_reader.indexed_fields
    ['name', 'gene_id']

    >>> # find all entries whose 'gene_id' matches 'nanos'
    >>> list(bb.search('gene_id', 'nanos'))
    [ list of matching SegmentChains/Transcripts ]



See also
--------
`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
    Descriptions of BED, GTF2, GFF3 and other text-based formats.
"""
import itertools
import struct
import sys
import warnings
import zlib
from collections import OrderedDict, Iterable
from plastid.genomics.roitools import GenomicSegment, SegmentChain, add_three_for_stop_codon
from plastid.readers.autosql import AutoSqlDeclaration
from plastid.util.io.binary import BinaryParserFactory, find_null_bytes
from plastid.util.unique_fifo import UniqueFIFO
from plastid.util.services.mini2to3 import ifilter, safe_bytes, safe_str
from plastid.util.services.decorators import skipdoc, deprecated
from plastid.util.services.exceptions import MalformedFileError, FileFormatWarning
from plastid.readers.autosql import AutoSqlDeclaration

from plastid.readers.bbifile cimport (
    bbiFile,
    bits32,
    bits64,
    lm,
    lmInit,
    lmCleanup,
    freeMem,
    _BBI_Reader,
    get_lm,
)

from plastid.genomics.roitools cimport GenomicSegment, SegmentChain
from plastid.genomics.c_common cimport (
    strand_to_str,
    str_to_strand,
    Strand,
    forward_strand,
    reverse_strand,
    unstranded,
    error_strand,
    _GeneratorWrapper,
)

from cpython.mem cimport PyMem_Malloc, PyMem_Free


#===============================================================================
# INDEX: BigBedReader
#===============================================================================

@skipdoc
class _FromBED_StrAdaptor(object):
    """Adaptor class to return strings from |BigBedReaders|.  Is internally
    called by BigBedReader when `return_type` is set to :class:`str`, because
    :class:`str` does not implement a ``from_bed`` method::

        >>> reader = BigBedReader(some_file, return_type=str)
    """
    @staticmethod
    def from_bed(inp):
        """Dummy method. Returns strings as themselves, instead of parsing a
        `BED`_ line

        Parameters
        ----------
        inp : str
            line of `BED`_ formatted text

        Returns
        -------
        str
            `inp`
        """
        return inp

cdef class BigBedReader(_BBI_Reader):
    """
    BigBedReader(filename, return_type = SegmentChain, add_three_for_stop = False, maxmem = 0)

    Reader for `BigBed`_ files. This class is useful for both iteration over
    genomic features one-by-one (like a reader), as well as random access to
    genomic features that overlap a region of interest (like a |GenomeHash|).

    Examples
    --------
    Iterate over all features in a `BigBed`_ file::

        >>> my_reader = BigBedReader("some_file.bb")
        >>> for feature in my_reader:
        >>>    pass # do something with each feature

    `BigBed`_ files can be accessed as dictionaries. To find features
    overlapping a region of interest::

        >>> roi = GenomicSegment("chrI", 0, 100000, "+")
        >>> for feature in my_reader[roi]:
        >>>     pass # do something with that feature

    Find features overlapping a genomic region of interest `roi`, on either
    strand::

        >>> for feature in my_reader.get(roi, stranded=False):
        >>>     pass # do something with that feature


    Parameters
    ----------
    filename : str
        Path to `BigBed`_ file

    return_type : |SegmentChain| or subclass, optional
        Type of feature to return from assembled subfeatures (Default:
        |SegmentChain|)

    add_three_for_stop : bool, optional
        Some annotation files exclude the stop codon from CDS annotations. If set to
        `True`, three nucleotides will be added to the threeprime end of each
        CDS annotation, **UNLESS** the annotated transcript contains explicit
        stop_codon feature. (Default: `False`)

    maxmem : float, optional
        Maximum desired memory footprint for C objects, in megabytes.
        May be temporarily exceeded if large queries are requested. Does not
        include memory footprint of Python objects. (Default: 0, no limit)


    Attributes
    ----------
    extension_fields : OrderedDict
        Dictionary mapping custom field names to their descriptions, if any
        custom fields are present

    extension_types : OrderedDict
        Dictionary mapping custom field names to objects that parse their types
        from strings

    filename : str
        Path to open BigBed file

    num_records : int
        Number of features in file

    num_chroms : int
        Number of contigs or chromosomes in file

    chroms : dict
        Dictionary mapping chromosome names to sizes

    return_type : class implementing a :py:meth:`from_bed` method, or str
        Class of object to return (Default: |SegmentChain|)
    """
    def __init__(self,
                 filename,
                 return_type=None,
                 add_three_for_stop=False,
                 **kwargs
                 ):
        """
        Parameters
        ----------
        filename : str
            Path to `BigBed`_ file

        return_type : |SegmentChain| or subclass, optional
            Type of feature to return from assembled subfeatures (Default:
            |SegmentChain|)

        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations.
            If set to `True`, three nucleotides will be added to the threeprime
            end of each CDS annotation, **UNLESS** the annotated transcript
            contains explicit stop_codon feature. (Default: `False`)

        maxmem : float
            Maximum desired memory footprint for C objects, in megabytes.  May
            be temporarily exceeded if large queries are requested.  Does not
            include memory footprint of Python objects.  (Default: 0, no limit)
        """
        cdef:
            str autosql

        self._bbifile = bigBedFileOpen(safe_bytes(filename))

        self.total_fields         = self._bbifile.fieldCount
        self.bed_fields           = self._bbifile.definedFieldCount
        self.num_extension_fields = self.total_fields - self.bed_fields
        self.extension_fields     = OrderedDict()
        self.extension_types      = OrderedDict()

        if self.num_extension_fields > 0:
            autosql = self._get_autosql_str()
            try:
                asql_parser            = AutoSqlDeclaration(autosql)
                self.extension_fields  = OrderedDict(
                    list(asql_parser.field_comments.items())[-self.num_extension_fields:]
                )

                for fieldname in self.extension_fields:
                    self.extension_types[fieldname] = asql_parser.field_formatters[fieldname]

            except AttributeError:
                warnings.warn(
                    "Could not find or could not parse autoSql declaration in BigBed file '%s': %s"
                    % (self.filename, autosql), FileFormatWarning
                )
                for x in range(self.num_extension_fields):
                    fieldname = "custom_%s" % x
                    self.extension_fields[fieldname] = "no description"
                    self.extension_types[ fieldname] = str

        if return_type == None:
            self.return_type = SegmentChain
        elif return_type == str:
            self.return_type = _FromBED_StrAdaptor
        else:
            self.return_type = return_type

        self.add_three_for_stop = add_three_for_stop

    property num_chroms:
        """Number of chromosomes in the `BigBed`_ file"""
        def __get__(self):
            return len(self.c_chroms())

    property num_records:
        """Number of features in file"""
        def __get__(self):
            return bigBedItemCount(self._bbifile)

    property bed_fields:
        """Number of standard `BED`_ format columns included in file"""
        def __get__(self):
            return self._bbifile.definedFieldCount

    property extension_fields:
        """Dictionary of names and types extra fields included in `BigWig`_/`BigBed`_ file"""
        def __get__(self):
            return self.extension_fields

    property return_type:
        """Return type of reader"""
        def __get__(self):
            if self.return_type == _FromBED_StrAdaptor:
                return str
            else:
                return self.return_type

    property indexed_fields:
        """Names of indexed fields in BigBed file. These are searchable by `self.search`"""
        def __get__(self):
            cdef:
                slName *orig_name       = bigBedListExtraIndexes(self._bbifile)
                slName *extension_names = orig_name
                list    names           = []

            if orig_name != NULL:
                names.append(safe_str(extension_names.name))
                while extension_names.next != NULL:
                    extension_names = extension_names.next
                    names.append(safe_str(extension_names.name))

                slFreeList(orig_name)

            return names

    def __str__(self):
        return "<%s records=%s chroms=%s>" % (self.__class__.__name__, self.num_records, self.num_chroms)

    def __repr__(self):
        return str(self)

    def _get_autosql_str(self):
        """Fetch `autoSql`_ field definition string from `BigBed`_ file, if present

        Returns
        -------
        str
            autoSql-formatted string, or "" if no autoSql definition present
        """
        cdef:
            char *  c_asql = bigBedAutoSqlText(self._bbifile)
            str     p_asql

        try:
            if c_asql is NULL:
                p_asql = ""
            else:
                p_asql = safe_str(c_asql)
        finally:
            freeMem(c_asql)

        return p_asql

    cdef list _bigbedinterval_to_bedtext(self, bigBedInterval * iv, Strand strand = unstranded):
        """Convert `bigBedIntervals` to lists of BED-formatted text entries.

        Parameters
        ----------
        iv : *bigBedInterval
            Pointer to first bigBedInterval

        strand : Strand, optional
            Strand that must be overlapped for features to be returned
            (Default: `unstranded`, return all entries)


        Returns
        -------
        list
            List of BED-formatted text entries
        """
        cdef:
            str bed_row, rest
            list ltmp     = []
            dict chromids = self._chromids

        if chromids is None:
            self._define_chroms()
            chromids = self._chromids

        while iv != NULL:
            bed_row = "\t".join("%s" % X for X in [chromids[iv.chromId], iv.start, iv.end])

            if self.total_fields > 3: # if iv.rest is populated
                rest = safe_str(iv.rest)

                # if strand info is present and we care about it
                if self.bed_fields >= 6:
                    items = rest.split("\t")
                    ivstrand = str_to_strand(items[2])

                    # no overlap - restart loop at next iv
                    if ivstrand & strand == 0:
                        iv = iv.next
                        continue

                bed_row = "%s\t%s" % (bed_row, rest)

            ltmp.append(bed_row)
            iv = iv.next

        return ltmp

    def search(self, field_name, *values):
        """Search indexed fields in the `BigBed`_ file for records matching
        `value` See `self.indexed_fields` for names of indexed fields and
        `self.extension_fields` for descriptions of extension fields.

        Parameters
        ----------
        field_name : str
            Name of field to search

        *values : one or more str
            Value(s) to match. If multiple are given, records matching any
            value will be returned.


        Yields
        ------
        object
            `self.return_type` of matching record in the `BigBed`_ file


        Examples
        --------
        Find all entries matching a given gene ID::

            # open file
            >>> bb = BigBedFile("some_file.bb")

            # which fields are searchable?
            >>> bb.indexed_fields
            ['name', 'gene_id']

            # find all entries whose 'gene_id' matches 'nanos'
            >>> bb.search('gene_id', 'nanos')
            [ list of matching segmentchains ]

            # find all entries whose 'gene_id' matches 'nanos' or 'oskar'
            >>> bb.search('gene_id', 'nanos', 'oskar')
            [ list of matching segmentchains ]


        Raises
        ------
        IndexError
            If field `field_name` is not indexed
        """
        cdef:
            int              fieldIdx
            int            * idx = &fieldIdx
            bptFile        * bpt
            bigBedInterval * iv
            lm             * buf = self._get_lm()
            list             ltmp
            object           outfunc   = self.return_type.from_bed
            list             etypes    = list(self.extension_types.items())

            str              stmp
            bytes            val
            char**           vals
            int              n, num_vals

        if field_name not in self.indexed_fields:
            raise KeyError("BigBed file '%s' has no index named '%s'" % (self.filename, field_name))
        else:
            bpt = bigBedOpenExtraIndex(self._bbifile, safe_bytes(field_name), idx)

        if len(values) == 1:
            val  = safe_bytes(values[0])
            iv   = bigBedNameQuery(self._bbifile, bpt, fieldIdx, val, buf)
        else:
            num_vals = len(values)
            vals     = <char**> PyMem_Malloc(num_vals*sizeof(char*))
            if not vals:
                raise MemoryError(
                    "BigBedReader.search(%s, %s): could not allocate memory"
                    % (field_name, ",".join(values))
                )

            for n, stmp in enumerate(values):
                val = safe_bytes(stmp)
                vals[n] = val

            iv  = bigBedMultiNameQuery(self._bbifile, bpt, fieldIdx, vals, num_vals, buf)
            PyMem_Free(vals)

        ltmp = self._bigbedinterval_to_bedtext(iv)
        bptFileDetach(&bpt)

        if self.add_three_for_stop == True:
            return _GeneratorWrapper(
                        (add_three_for_stop_codon(outfunc(X, extra_columns=etypes)) for X in ltmp),
                        "BigBed entries"
                    )
        else:
            return _GeneratorWrapper(
                (outfunc(X, extra_columns=etypes) for X in ltmp),
                "BigBed entries"
            )

    def get(self, roi, bint stranded=True, bint check_unique=True):
        """Iterate over features that share genomic positions with a region of
        interest

        Note
        ----
        ``reader.get(roi)`` is an alternative syntax to ``reader[roi]``. It's
        only useful if setting `stranded` or `check_unique` to `False`.


        Parameters
        ----------
        roi : |SegmentChain| or |GenomicSegment|
            Query feature representing region of interest

        stranded : bool, optional
            If `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands. (Default: `True`)

        check_unique: bool, optional
            if `True`, assure that all results in generator are unique.
            (Default: `True`)


        Yields
        ------
        object
            `self.return_type` of each record in the `BigBed`_ file


        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef:
            SegmentChain chain

        if isinstance(roi, SegmentChain):
            chain = roi
        elif isinstance(roi, GenomicSegment):
            chain = SegmentChain(roi)
        else:
            raise TypeError("BigBedReader.get(): Query interval must be a GenomicSegment or SegmentChain")

        return self._c_get(chain, stranded, check_unique=check_unique)

    # TODO: direct C/Cython route to SegmentChain.from_bed
    # NB- no cache layer, which we  had in pure Python implementation (below)
    # will this be fast enough for repeated queries over the same region?
    # need to test
    cdef _GeneratorWrapper _c_get(
        self,
        SegmentChain chain,
        bint stranded=True,
        bint check_unique=True,
        lm *my_lm = NULL
        ):
        """c-layer implementation of :meth:`BigBedReader.get`

        Parameters
        ----------
        roi : |SegmentChain| or |GenomicSegment|
            Query feature representing region of interest

        stranded : bool, optional
            If `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands. (Default: `True`)

        check_unique : bool, optional
            If `True`, assure that all results in generator are unique.
            (Default: `True`)

        my_lm : lm, optional
            If not NULL, use this pool of local memory instead of the
            |BigBedReader|'s.  (Default: NULL)


        Yields
        ------
        object
            `self.return_type` of each record in the `BigBed`_ file
        """
        cdef:
            bigBedInterval * iv
            list             ltmp      = []
            GenomicSegment   span      = chain.spanning_segment
            GenomicSegment   roi
            str              chrom     = span.chrom
            Strand           strand    = unstranded
            lm*              buf #       = self._get_lm()
            Strand           ivstrand
            object           outfunc   = self.return_type.from_bed
            list             etypes    = list(self.extension_types.items())

        if my_lm != NULL:
            buf = my_lm
        else:
            buf = self._get_lm()

        if stranded is True:
            strand = span.c_strand

        for roi in chain:
            iv = bigBedIntervalQuery(self._bbifile,
                                     safe_bytes(span.chrom),
                                     roi.start,
                                     roi.end,
                                     0,
                                     buf)
            ltmp.extend(self._bigbedinterval_to_bedtext(iv, strand=strand))

        # filter for uniqueness
        if check_unique == True:
            ltmp = sorted(set(ltmp))

        if self.add_three_for_stop == True:
            return _GeneratorWrapper((
                    add_three_for_stop_codon(outfunc(X, extra_columns=etypes)) for X in ltmp
                    ),
                "BigBed entries"
                )
        else:
            return _GeneratorWrapper((outfunc(X, extra_columns=etypes) for X in ltmp), "BigBed entries")

    def __getitem__(self, roi):
        """Iterate over features that share genomic positions with a region of
        interest, on same strand.  Unstranded features, if present, are
        considered to overlap both `rois` on any strand.


        Parameters
        ----------
        roi : |SegmentChain| or |GenomicSegment|
            Query feature representing region of interest

        Yields
        ------
        object
            `self.return_type` of each record in the `BigBed`_ file


        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        return self.get(roi, stranded=True)

    def __iter__(self):
        """Iterate over all features in `BigBed`_ file

        Yields
        ------
        object
            Object of ``self.return_type``, |SegmentChain| or one of its
            subclasses
        """
        return _GeneratorWrapper(BigBedIterator(self, maxmem=self._maxmem), "BigBed records")


# can't be cdef'ed or cpdef'ed due to yield
#
# This would probably be faster if we iterated through the cirTree leaf nodes
# directly, since they overlap chromosome boundaries. But, as long as the
# number of chromosomes/contigs in the file is low, this shouldn't be too big
# a deal.
#
# But, cirTree.h, crTree.h, and BigBed.h don't give a convenient way to get
# all the data from a node, so doing so would result in lots of reimplementation.
def BigBedIterator(BigBedReader reader, maxmem=0):
    """BigBedIterator(reader, maxmem = 0)

    Iterate over records in the `BigBed`_ file, sorted lexically by chromosome
    and position.

    Parameters
    ----------
    reader : |BigBedReader|
        Reader to iterate over

    maxmem : float
        Maximum desired memory footprint for C objects, in megabytes.  May be
        temporarily exceeded if large queries are requested.  Does not include
        memory footprint of Python objects.  (Default: 0, no limit)

    Yields
    ------
    object
        reader.return_type of `BED`_ record

    Raises
    ------
    MemoryError
        If memory cannot be allocated
    """
    cdef:
        list           chromsizes  = sorted(reader.c_chroms().items(), key = lambda x: x[0].lower())
        long           chromlength
        str            chrom
        lm             *buf = lmInit(0)

    if not buf:
        raise MemoryError("BigBedIterator: could not allocate local memory")

    for chrom, chromlength in chromsizes:
        query = SegmentChain(GenomicSegment(chrom, 0, chromlength, "."))
        buf = get_lm(my_lm=buf, maxmem=maxmem)
        if not buf:
            raise MemoryError("BigBedIterator: could not allocate local memory")

        for n, roi in enumerate(reader._c_get(query, stranded=False, check_unique=False, my_lm=buf)):
            yield roi

    lmCleanup(&buf)


#===============================================================================
# INDEX: Factories for various record formats
#===============================================================================

HeaderFactory = BinaryParserFactory(
    "BigBedHeader",
    "IHH3QHHQQIQ",
    [
        "magic",
        "version",
        "zoom_levels",
        "bplus_tree_offset",
        "full_data_offset",
        "r_tree_offset",
        "field_count",
        "bed_field_count",
        "autosql_offset",
        "total_summary_offset",
        "uncompressed_buffer_size",
        "reserved",
    ],
)
"""Parse headers for BigWig and `BigBed`_ files"""

ZoomHeaderFactory = BinaryParserFactory(
    "ZoomHeader",
    "2H2Q",
    [
        "reduction_level",
        "reserved",
        "data_offset",
        "index_offset",
    ],
)
"""Parse Zoomlevel headers in BigWig files"""


TotalSummaryFactory = BinaryParserFactory(
    "TotalSummary",
    "5Q",
    [
        "bases_covered",
        "min_val",
        "max_val",
        "sum_data",
        "sum_squares",
    ],
)
"""Parse 'Total Summary' tables in BigWig and `BigBed`_ files"""

BPlusTreeHeaderFactory = BinaryParserFactory(
    "BPlusTreeHeader",
    "4I2Q",
    [
        "magic",
        "block_size",
        "key_size",
        "val_size",
        "num_chroms",
        "reserved",
    ],
)
"""Parse headers of |BPlusTree|"""

BPlusTreeNodeFormatFactory  = BinaryParserFactory(
    "BPlusTreeNodeFormat",
    "?BH",
    ["is_leaf", "reserved", "child_count"],
)
"""Determine format of BPlus Tree nodes in a datablock of a `BigBed`_ or BigWig file"""

RTreeHeaderFactory = BinaryParserFactory(
    "RTreeHeader",
    "IIQ4IQII",
    [
        "magic",
        "block_size",
        "num_nodes",
        "start_chrom_id",
        "start_base",
        "end_chrom_id",
        "end_base",
        "end_file_offset",
        "items_per_leaf",
        "reserved",
    ]
)
"""Parse headers of |RTree|"""

RTreeNodeFormatFactory = BinaryParserFactory(
    "RTreeNodeFormat",
    "BBH",
    ["is_leaf", "reserved", "count"],
    )
"""Determine format of R Tree nodes in a data block of a `BigBed`_ or BigWig file"""


RTreeLeafFactory = BinaryParserFactory(
    "RTreeLeaf",
    "4I2Q",
    [
        "start_chrom_id",
        "start_base",
        "end_chrom_id",
        "end_base",
        "data_offset",
        "data_size",
    ],
)
"""Parse leaves of |RTree|"""

RTreeNonLeafFactory = BinaryParserFactory(
    "RTreeNonLeaf",
    "4IQ",
    [
        "start_chrom_id",
        "start_base",
        "end_chrom_id",
        "end_base",
        "child_data_offset",
    ],
)
"""Parse non-leaf nodes of |RTree|"""
