#!/usr/bin/env python
"""|BigBedReader|, a parser for `BigBed`_ files. In contrast to `BED`_, `GTF2`_,
and `GFF3`_ files, `BigBed`_ files are binary, indexed, and randomly-accessible.
This has several interesting implications:

    - |BigBedReader| can be used to iterate over records, like a reader, **or**
      to fetch records that cover a region of interest, in the manner of a |GenomeHash|
    
    - Because `BigBed`_ files are randomly accessible, their records don't need
      to be loaded into memory to be parsed, which results in substantial memory
      savings for large genomes. This, however, comes at a cost in speed, 
      as `BigBed`_ files are highly compressed.


Examples
--------
Iterate over all features in a `BigBed`_ file::

    >>> my_reader = BigBedReader("some_file.bb")
    >>> for feature in my_reader:
            pass # do something with each feature


Instead, find features overlapping a specific feature or region of interest `roi`,
on the same strand as the `roi`::

    >>> roi = GenomicSegment("chrI",0,100000,"+")
    >>> for feature in my_reader[roi]:
            pass # do something with that feature
            ...

    
Find features overlapping a genomic region of interest `roi`,
on either strand::

    >>> for feature in my_reader.__getitem__(roi,stranded=False):
            pass # do something with that feature



See also
--------
|BigBedReader|
    Class documentation for |BigBedReader|

`Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_
    Description of BigBed and BigWig formats. Especially see supplemental data.

`UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
    Descriptions of BED, GTF2, GFF3 and other text-based formats.
"""
import struct
import zlib
import itertools
import sys
import warnings
from collections import OrderedDict
from plastid.genomics.roitools import GenomicSegment, SegmentChain
from plastid.readers.common import add_three_for_stop_codon
from plastid.readers.autosql import AutoSqlDeclaration
from plastid.util.io.binary import BinaryParserFactory, find_null_bytes
from plastid.util.io.openers import NullWriter
from plastid.util.unique_fifo import UniqueFIFO
from plastid.util.services.mini2to3 import ifilter
from plastid.util.services.decorators import skipdoc
from plastid.util.services.exceptions import MalformedFileError, FileFormatWarning

#===============================================================================
# INDEX: BigBedReader
#===============================================================================

@skipdoc
class _FromBED_StrAdaptor(object):
    """Adaptor class to return strings from |BigBedReaders|.
    Is internally called by BigBedReader when `return_type` is set to :class:`str`,
    because :class:`str` does not implement a ``from_bed`` method::

        >>> reader = BigBedReader(some_file,return_type=str)
    """
    @staticmethod
    def from_bed(inp):
        """Dummy method. Returns strings as themselves, instead of parsing a `BED`_ line

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


class BigBedReader(object):
    """Reader for `BigBed`_ files. This class is useful for both iteration
    over genomic features one-by-one (like a reader), as well as random access to
    genomic features that overlap a region of interest (like a |GenomeHash|).  
    
    
    Examples
    --------
    Iterate over all features in a BigBed file::
    
        >>> my_reader = BigBedReader("some_file.bb")
        >>> for feature in my_reader:
                pass # do something with each feature


    Instead, find features overlapping a specific region of interest `roi`,
    on the same strand as the `roi`, and save these as a list::
    
        >>> roi = GenomicSegment("chrI",0,100000,"+")
        >>> overlapping_features = list(my_reader[roi])
    

    Find features overlapping a genomic region of interest `roi`,
    on either strand::
    
        >>> for feature in my_reader.__getitem__(roi,stranded=False):
                pass # do something with that feature
                ...
    
    
    Attributes
    ----------
    custom_fields : OrderedDict
        Dictionary mapping custom field names to their descriptions,
        if any custom fields are present

    header : dict
        Dictionary containing header information

    filename : str
        Path to open BigBed file

    num_records : int
        Number of features in file
    
    num_chroms : int
        Number of contigs or chromosomes in file
    
    chrom_sizes : dict
        Dictionary mapping chromosome names to sizes
    
    return_type : class implementing a :py:meth:`from_bed` method, or str
        Class of object to return (Default: |SegmentChain|)


    Notes
    -----
    See `Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_ for 
    detailed description of the structures of the BigBed format, B+ Tree,
    and R trees. 
    """
    def __init__(self,
                 filename,
                 base_record_format="III",
                 return_type=SegmentChain,
                 memorize_r_tree=False,
                 add_three_for_stop=False,
                 cache_depth=5,
                 printer=NullWriter()
                 ):
        """Create a BigBedReader
        
        Parameters
        ----------
        filename : str
            String indicating path to `BigBed`_ file (*not* open filehandle)
        
        base_record_format : str, optional
            Format string for :py:func:`struct.unpack`, excluding endian-ness prefix
            and any notion of a null-terminated string (Default: "III")
        
        return_type : str or class implementing a :py:meth:`from_bed` method, optional
            Class of object to return (Default: |SegmentChain|)
        
        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations. If set to
            `True`, and the `BigBed`_ file annotates transcripts, three nucleotides
            will be added to the threeprime end of each CDS annotation.
            (Default: `False`)
                        
        cache_depth : int, optional
            Number of previously-fetched data blocks to keep in memory.
            Decrease this number to reduce memory usage. Increase it to speed up
            repeated fetches to nearby genomic regions. (Default: `5`)

        memorize_r_tree : bool, optional
            If *True*, cache entire |RTree| index into memory for faster lookup
            (faster for small files, use big memory for large files. Default: `False`)
            
        printer : file-like, optional
            Filehandle or sys.stderr-like for logging (Default: |NullWriter|)            
        """
        self.filename = filename
        self.fh = open(filename,"rb")
        if return_type == str:
            self.return_type = _FromBED_StrAdaptor 
        else:
            self.return_type = return_type

        self.base_record_format = base_record_format
        self.add_three_for_stop = add_three_for_stop

        # whether or not file is byteswapped
        # this is re-detected below when the file in _parse_header
        # default: little-endian
        self._byte_order = "<"
        
        self.header      = self._parse_header()
        self.num_records = self._count_records()
        
        self._num_custom_fields  = self.header["field_count"] - self.header["bed_field_count"]
        self.custom_fields = {}
        
        if self._num_custom_fields > 0:
            autosql = self._get_autosql_str()
            try:
                self.autosql_parser = AutoSqlDeclaration(autosql)
                self.custom_fields  = OrderedDict(list(self.autosql_parser.field_comments.items())[-self._num_custom_fields:])
            except AttributeError:
                warnings.warn("Could not find or could not parse autoSql declaration in BigBed file '%s': %s" % (self.filename,autosql),FileFormatWarning)
                self.custom_fields  = OrderedDict([("custom_%s" % X,"no description") for X in range(self._num_custom_fields)])
                self.autosql_parser = lambda x: OrderedDict(zip(self.custom_fields,x.split("\t")[-self._num_custom_fields:]))
        
        self.bplus_tree = BPlusTree(self.filename,self.header["bplus_tree_offset"],self._byte_order)
        self.bplus_tree_header  = self.bplus_tree.header
        self.num_chroms  = self.bplus_tree.num_chroms
        self.chrom_sizes = self.bplus_tree.chrom_sizes

        self.r_tree = RTree(self.filename,
                            self.bplus_tree,
                            self.header["r_tree_offset"],
                            self._byte_order,
                            memorize=memorize_r_tree)
        
        self.fifo = UniqueFIFO(cache_depth)
        self.fifo_dict = {}

    def close(self):
        """Close all open pointers to `BigBed`_ file"""
        self.fh.close()
        self.r_tree.fh.close()
        self.bplus_tree.fh.close()
        
    def __del__(self):
        self.close()

    def __str__(self):
        return "<%s records=%s chroms=%s>" % (self.__class__.__name__,self.num_records,self.num_chroms)

    def __repr__(self):
        return str(self)

    def _get_autosql_str(self):
        """Fetch `autoSql`_ field definition string from `BigBed`_ file, if present
        
        Returns
        -------
        str
            autoSql-formatted string
        """
        if self.header["autosql_offset"] > 0:
            self.fh.seek(self.header["autosql_offset"])
            autosql_len = self.header["total_summary_offset"] - self.header["autosql_offset"]
            autosql_fmt = "%s%ss" % (self._byte_order,autosql_len)
            autosql_str, = struct.unpack(autosql_fmt,self.fh.read(autosql_len))
            return str(autosql_str.decode("ascii")).strip("\0")
        else:
            return ""
            
    def _count_records(self):
        """Counts number of features in `BigBed`_ file
        
        Returns
        -------
        int
            Number of features in `BigBed`_ file
        """
        self.fh.seek(self.header["full_data_offset"])
        return struct.unpack(self._byte_order+"Q",self.fh.read(8))[0]

    def _parse_header(self):
        """Parse first 64 bytes of `BigBed`_ file, and determine indices
        of file metadata. 
        
        Header table information from :cite:`Kent2010`, Supplemental table 5:
        
        =========================  ==== ====  =================================================
        Field                      Size Type   Summary
        =========================  ==== ====  =================================================
        magic                      4    uint   0x8789F2EB
        version                    2    uint   File version 3?
        zoom_levels                2    uint   Number of zoom summary resolutions
        bplus_tree_offset          8    uint   Offset to chr B+ tree index
        full_data_offset           8    uint   Offset to main data. dataCount
        r_tree_offset              8    uint   Offset to R tree index of items
        field_count                2    uint   Number of fields in BED file
        bed_field_count            2    uint   Number of fields that are pre-defined BED fields
        autosql_offset             8    uint   Offset to zero-terminated string with .as spec. 0 if no autoSql string present
        total_summary_offset       8    uint   Offset to overall file summary data block
        uncompressed_buffer_size   4    uint   Maximum size decompression buffer needed (files v3 and later)
        reserved                   8    uint   Reserved for future expansion
        =========================  ==== ====  =================================================
        
        Returns
        -------
        dict
            Dictionary mapping header field names to values 
        """
        self.fh.seek(0)
        items = HeaderFactory(self.fh,self._byte_order)

        if items["magic"] != 0x8789F2EB:
            self._byte_order = ">"
            self.fh.seek(0)
            items = HeaderFactory(self.fh,self._byte_order)
        
        if items["magic"] != 0x8789F2EB:
            raise MalformedFileError(self.filename,"Could not determine byte order of BigBed file. Expected magic number to be '%x', got '%x'." % (0x8789F2EB,self.header["magic"]))
        
        return items
    
    def _iterate_over_chunk(self,data_offset,data_size,null=b"\x00"):
        """Iterate over records in a portion of the `BigBed`_ file
        
        Parameters
        ----------
        data_offset : int
            Beginning of compressed record block in `BigBed`_ file
        
        data_size : int
            End of compressed record block in `BigBed`_ file
        
        null : str, optional
            Null character. Default: `'\\x00'`
        
        Yields
        ------
        object
            Object of `self.return_type`, usually |SegmentChain| or one of its subclasses
        """
        base_size = struct.calcsize(self._byte_order+self.base_record_format)
        
        if data_offset in self.fifo:
            self.fifo.append(data_offset) # this line looks funny here, but just moves data_offset to end of FIFO
            for my_obj in self.fifo_dict[data_offset]:
                yield my_obj
        else:
            self.fifo.append(data_offset) # here, it actually adds data_offset to end of FIFO
            
            # remove any data that should be kicked out of FIFO
            # this would be in the keyset to self.fifo_dict ubt not in self.fifo
            for k in set(self.fifo_dict.keys()) - set(self.fifo._elements):
                self.fifo_dict.pop(k)

            self.fifo_dict[data_offset] = []
        
            self.fh.seek(data_offset)
            raw_data = zlib.decompress(self.fh.read(data_size))
            if self.header["field_count"] > 3:
                null_indices = find_null_bytes(raw_data,null=null)
                last_index = 0
                while (null_indices > last_index + base_size).sum() > 0:
                    # find first index that doesn't overlap numerical data
                    end_index = null_indices[(null_indices > last_index + base_size).argmax()]
                
                    # get bytes covering next record
                    bytestr = raw_data[last_index:end_index]
                    str_bytes = end_index - last_index - base_size
                    new_format = "%s%s%ss" % (self._byte_order,self.base_record_format,str_bytes)
                    
                    # reset starting index
                    last_index = end_index + 1
                    
                    # parse record
                    chrom_id, chrom_start, chrom_end, remaining_cols = struct.unpack(new_format,bytestr)
                    
                    # 2.x returns str, 3.x returns bytes explicit cast here.
                    if isinstance(remaining_cols,bytes):
                        remaining_cols = remaining_cols.decode("ascii")
                    if sys.version_info < (3,) and isinstance(remaining_cols,unicode):
                        remaining_cols = str(remaining_cols.decode("ascii"))

                    chrom_name = self.bplus_tree.chrom_id_name[chrom_id]
                    items      = remaining_cols.split("\t")
                    bed_items  = items[:self.header["bed_field_count"] - 3]
                    bed_line   = "\t".join([chrom_name,str(chrom_start),str(chrom_end)] + bed_items)
                    return_obj = self.return_type.from_bed(bed_line)
        
                    # if file contains custom fields, parse these and put them in attr dict
                    whole_bedplus_line = "\t".join([chrom_name,str(chrom_start),str(chrom_end)])+"\t"+remaining_cols
                    if self._num_custom_fields > 0:
                        return_obj.attr.update(list(self.autosql_parser(whole_bedplus_line).items())[-self._num_custom_fields:])
                    
                    if self.add_three_for_stop == True:
                        return_obj = add_three_for_stop_codon(return_obj)
                    
                    self.fifo_dict[data_offset].append(return_obj)
                    yield return_obj
            else:
                last_index = 0
                fmt_str = "%s%ss" % (self._byte_order,self.base_record_format)
                calcsize = struct.calcsize(fmt_str) 
                while last_index < len(raw_data):
                    bytestr = raw_data[last_index:last_index+calcsize]
                    chrom_id, chrom_start, chrom_end, _ = struct.unpack(fmt_str,bytestr)
                    chrom_name = self.bplus_tree.chrom_id_name[chrom_id]
                    bed_line   = "\t".join([chrom_name,str(chrom_start),str(chrom_end)])
                    return_obj = self.return_type.from_bed(bed_line)
                    last_index = last_index+calcsize
                    
                    # don't need to add three for stop when field count < 3,
                    # because no stop codon
                    yield return_obj
            
    def __getitem__(self,roi,stranded=True):
        """Iterate over features that overlap a region of interest
        
        Parameters
        ----------
        roi : |SegmentChain| or |GenomicSegment|
            Query feature representing region of interest
        
        stranded : bool, optional
            if `True`, retrieve only features on same strand as query feature.
            Otherwise, retrieve features on both strands. (Default: `True`)
            
        
        Yields
        ------
        object
            Object of `self.return_type`, |SegmentChain| or one of its subclasses
        
        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        if isinstance(roi,SegmentChain):
            segchain = roi
        elif isinstance(roi,GenomicSegment):
            segchain = SegmentChain(roi)
        else:
            raise TypeError("Query interval must be a GenomicSegment or SegmentChain")
        
        if stranded:
            overlap_fn = SegmentChain.overlaps
        else:
            overlap_fn = SegmentChain.unstranded_overlaps
            
        return ifilter(lambda x: overlap_fn(segchain,x),
                                 itertools.chain.from_iterable((self._iterate_over_chunk(file_offset,byte_length) \
                                                                for (file_offset,byte_length) \
                                                                in self.r_tree[segchain.spanning_segment])))
    
    def __iter__(self):
        """Generator that iterates over all features in `BigBed`_ file
        
        Yields
        ------
        object
            Object of ``self.return_type``, |SegmentChain| or one of its subclasses
        """
        return itertools.chain.from_iterable((self._iterate_over_chunk(*NODE) for NODE in self.r_tree))

#===============================================================================
# INDEX: BPlusTree parser
#===============================================================================

class BPlusTree(object):
    """Decode B+ Trees, which are used to describe chromosomes and contigs
    in `BigBed`_ and BigWig files.

    See `Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_ for 
    detailed description of the structures of the `BigBed`_ format, B+ Tree,
    and R trees.
    
    
    Attributes
    ----------
    filename : str
        Path to `BigBed`_ file

    header : dict
        Dictionary describing |BPlusTree|
    
    header_offset : int
        Offset in BigBed file to BPlus Tree header data
    
    tree_offset : int
        Offset in BigBed file to root node info block of BPlus tree
        
    num_chroms : int
        Number of chromosomes described y the BPlus tree

    chrom_sizes : dict
        Dictionary mapping chromosome names to sizes in base pairs
    
    chrom_id_name : dict
        Dictionary mapping internal chromosome IDs to chromosome names
    
    chrom_name_id : dict
        Dictionary mapping chromosome names to internal chromosome IDs
    
    LeafInfoFactory : |BinaryParserFactory|
        Binary parser for leaf nodes
    
    NonLeafInfoFactory : |BinaryParserFactory| 
        Binary parser for non-leaf nodes
    """
    
    def __init__(self,filename,start_offset,byte_order="<"):
        """Create a |BPlusTree|
        
        Parameters
        ----------
        filename : str
            String indicating path to `BigBed`_ file (*not* open filehandle)
        
        start_offset : int
            Offset, in bytes, to BPlus Tree Header
        
        byte_order : str
            Character indicating endian-ness of data (default: `'<'` for little-endian)
        """
        self._byte_order = byte_order
        self.filename = filename
        self.fh = open(filename,"rb")
        self.header_offset = start_offset
        self.header      = self._parse_header()

        if self.header["magic"] != 0x78CA8C91:
            raise MalformedFileError(self.filename,"Could not determine byte order of B+ Tree. Expected magic number to be %x, got %x." % (0x78CA8C91,self.header["magic"]))

        self.tree_offset = start_offset + BPlusTreeHeaderFactory.calcsize()
        self.num_chroms = self.header["num_chroms"]
        self.LeafInfoFactory = BinaryParserFactory("BTreeLeafInfo",
                                                   "%ssII" % self.header["key_size"],
                                                  ["chrom_name",
                                                   "chrom_id",
                                                   "chrom_size"])
        self.NonLeafInfoFactory = BinaryParserFactory("BTreeNonLeafInfo",
                                                      "%ssQ" % self.header["key_size"],
                                                     ["chrom_name",
                                                      "child_offset"])
        
        self.chrom_sizes   = OrderedDict()
        self.chrom_id_name = OrderedDict()
        self.chrom_name_id = OrderedDict()

        self._get_chrom_data()

    def __del__(self):
        self.fh.close()
    
    def _parse_header(self):
        """Parses BPlus Tree Header in `BigBed`_ file
        
        Header table information from :cite:`Kent2010`, Supplemental table 8:
        
        =====================  =====  =====  ========================================
        Field                  Size   Type   Summary
        =====================  =====  =====  ========================================
        magic                  4      uint   0x78CA8C91
        block_size             4      uint   Number of children per block
        key_size               4      uint   Number of significant bytes per key,
                                             the minimum number of bytes needed
                                             to distinguish chromosome names
        val_size               4      uint   Size of value being indexed. Currently 8
        num_contigs            8      uint   Number of chromosomes or contigs
        reserved               8      uint   Reserved for future expansion        
        =====================  =====  =====  ========================================
        
        Returns
        -------
        dict
            Dictionary containing header info
        """
        self.fh.seek(self.header_offset)
        return BPlusTreeHeaderFactory(self.fh,self._byte_order)
    
    def _get_chrom_data(self):
        """Retrieves chromosome sizes and names from BPlus tree, and hashes
        these as dictionaries keyed on chromosome name or ID.
        
        Populates the following:
        
        `self.chrom_id_name`
            Maps chromosome IDs to names
        
        `self.chrom_name_id`
            Maps chromosome names to IDs
        
        `self.chrom_sizes`
            Maps chromosome names to sizes in basepairs
        """
        for chrom_id,chrom_name,chrom_size in self._walk_tree():
            self.chrom_id_name[chrom_id]   = chrom_name
            self.chrom_name_id[chrom_name] = chrom_id
            self.chrom_sizes[chrom_name]   = chrom_size
        
    def _walk_tree(self,start_offset=None):
        """Exhaustively traverses BPlus tree, starting at the node starting at `start_offset`
        
        Parameters
        ----------
        start_offset : int
            Offset of node block in file
        
        Returns
        -------
        list
            List of tuples of `(chrom_id,chrom_name,chrom_size)`, in order
            of traversal depth-first, left-to-right
        """
        if start_offset is None:
            start_offset = self.tree_offset
            
        fh = open(self.filename,"rb")
        fh.seek(start_offset)
        
        chrom_info = []
        
        node_info = BPlusTreeNodeFormatFactory(fh)
        if node_info["is_leaf"] == True:
            for _ in range(node_info["child_count"]):
                new_node = self.LeafInfoFactory(fh)
                chrom_name = new_node["chrom_name"].strip("\x00")
                chrom_info.append((new_node["chrom_id"],chrom_name,new_node["chrom_size"]))
        else:
            for _ in range(len(node_info["child_count"])):
                new_node = self.NonLeafInfoFactory(fh)
                chrom_info.extend(self.get_chrom_data2(start_offset=new_node["child_offset"]))
        
        fh.close()
        return chrom_info
        
#===============================================================================
# INDEX: R tree parser
#===============================================================================

class RTree(object):
    """Decode R Trees, which index genomic coordinates to file positions in `BigBed`_
    and BigWig files

    See `Kent2010 <http://dx.doi.org/10.1093/bioinformatics/btq351>`_ for 
    detailed description of the structures of the BigBed format, B+ Tree,
    and R trees. 

    Attributes
    ----------
    filename : str
        Path to `BigBed`_ file
    
    header : dict
        Dictionary describing |RTree|

    leaf_data : dict
        Dictionary mapping record offsets of leaf nodes to tuples of their
        `(data_offset,data_size)`. This dictionary is populated via lazy
        evaluation
        
    header_offset : int
        Offset in `BigBed`_ file to |RTree| Tree header data
    
    tree_offset : int
        Offset in `BigBed`_ file to root node info block of |RTree| tree

    bplus_tree : |BPlusTree|
        |BPlusTree| associated with the `BigBed`_ file
    """
    
    def __init__(self,
                 filename,
                 bplus_tree,
                 start_offset,
                 byte_order,
                 memorize=False):
        """Create an R Tree
        
        Parameters
        ----------
        filename : str
            String indicating path to `BigBed`_ file (*not* open filehandle)
        
        start_offset : int
            Offset, in bytes, to |RTree| header
        
        byte_order : str
            Character indicating endian-ness of data (default: `'<'` for little-endian)
            
        memorize : bool
            If `True`, cache entire tree topology into memory for faster lookup
            (faster for small files, uses big memory for large files. Default: `False`)
        """     
        self.filename      = filename
        self.bplus_tree    = bplus_tree
        self._byte_order   = byte_order
        self.memorize      = memorize
        self.header_offset = start_offset
        self.tree_offset   = start_offset + RTreeHeaderFactory.calcsize()

        self.fh            = open(filename,"rb")
        self.header        = self._parse_header()
        if self.header["magic"] != 0x2468ACE0:
            raise MalformedFileError(self.filename,"Could not determine byte order of R tree. Expected magic number to be %x. Got %x." % (0x2468ACE0,self.header["magic"]))
        
        # to be populated as needed
        self.leaf_data              = {}
        self.node_parent_offsets    = {}
        self.node_child_offsets     = {}
        self.node_genome_boundaries = {}

    def __del__(self):
        self.fh.close()

    def _parse_header(self):
        """Parses |RTree| Header in `BigBed`_ file
        
        Header table information from :cite:`Kent2010`, Supplemental table 8:
        
        ===================  ======  ======  =================================================
        Field                Size    Type    Summary
        ===================  ======  ======  =================================================
        magic                  4      uint   0x2468ACE0
        block_size             4      uint   Number of children per block
        num_nodes              8      uint   Number of [leaf?] nodes
        start_chrom_id         4      uint   ID of first chromosome in index
        start_base             4      uint   Position of first base in index
        end_chrom_id           4      uint   ID of last chromosome in index
        end_base               4      uint   Position of last base in index
        end_file_offset        8      uint   Position in file where data begin indexed ends
        items_per_leaf         4      uint   Number of items pointed to by leaves of index
        reserved               4      uint   Reserved
        ===================  ======  ======  =================================================
        
        Returns
        -------
        dict
            Dictionary containing header info
        """
        self.fh.seek(self.header_offset)
        return RTreeHeaderFactory(self.fh,self._byte_order)
    
    def __getitem__(self,roi,start_node_offset=None):
        """Search |RTree| and return addresses of data blocks covering a region of interest

        Parameters
        ----------
        roi : |SegmentChain| or |GenomicSegment|
            Query feature representing region of interest

        start_node_offset : int, optional
            Start offset of first node block (Default: top of |RTree|) 
                     
        Returns
        -------
        list
            List of tuples of `(file_offset, bytes_to_read)`  in `BigBed` file
            specifying data blocks of find records that should be checked for overlap
            with `roi`
         
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        if isinstance(roi,SegmentChain):
            iv = roi.iv
        elif isinstance(roi,GenomicSegment):
            iv = roi
        else:
            raise TypeError("Query interval must be a GenomicSegment or SegmentChain")

        if start_node_offset is None:
            start_node_offset = self.tree_offset
                     
        fh = open(self.filename,"rb")
        fh.seek(start_node_offset)
        node_info = RTreeNodeFormatFactory(fh,self._byte_order)
        
        chrom_id = self.bplus_tree.chrom_name_id.get(iv.chrom,None)
        if chrom_id is None:
            return []

        ltmp = []
        
        for _ in range(node_info["count"]):
            
            if node_info["is_leaf"] == True:
                new_node = RTreeLeafFactory(fh)
                # if leaf overlaps, append (data_offset, data_size) to list
                if self.node_overlaps_roi(new_node,chrom_id,iv.start,iv.end):
                    ltmp.append((new_node["data_offset"],new_node["data_size"]))
                
            else:
                new_node = RTreeNonLeafFactory(fh)
                if self.node_overlaps_roi(new_node,chrom_id,iv.start,iv.end):
                    ltmp.extend(self.__getitem__(iv,start_node_offset=new_node["child_data_offset"]))

        fh.close()
        return ltmp

    @staticmethod
    def node_overlaps_roi(node,roi_chrom_id,roi_start_base,roi_end_base):
        """Determines whether or not an |RTree| node overlaps a region of interest (ROI)
        
        Parameters
        ----------
        node : RTreeLeaf or RTreeNonLeaf
            Query node
        
        roi_chrom_id: int
            Integer corresponding to chromosome ID for ROI
        
        roi_start_base : int
            Coordinate of leftmost genomic position of ROI
        
        roi_end_base : int
            Coordinate of rightmost genomic position of ROI
        
        Returns
        -------
        bool
            `True` if `node` overlaps the ROI. `False` otherwise
        """
        if node["start_chrom_id"] > roi_chrom_id:
            return False
        elif node["start_chrom_id"] == roi_chrom_id:
            if node["start_base"] >= roi_end_base:
                return False
            else:
                # overlap if end base is after start; or if end base is on successive chromosome
                if node["end_base"] > roi_start_base or node["end_chrom_id"] > roi_chrom_id:
                    return True
                else:
                    return False
        else: #start_chrom_id < roi_chrom_id
            if node["end_chrom_id"] < roi_chrom_id:
                return False
            elif node["end_chrom_id"] == roi_chrom_id:
                if node["end_base"] <= roi_start_base: # chromosome ends before region starts
                    return False
                else:
                    # last chromosome in block ends after region starts
                    # first chromosome begins earlier
                    return True
            else:
                return True
    
    def __iter__(self):
        """Return an iterator over the leaf nodes of the |RTree|"""
        # find leaves if we haven't already cached them
        if len(self.leaf_data) == 0:
            self._find_leaves()
        
        return (items for items in sorted(self.leaf_data.values()))
        
    def _find_leaves(self,start_node_offset=None):
        """Search |RTree| exhaustively and return addresses of data blocks of all
        leaves. Leaf data is stored in:
        
            `self.leaf_data`
                Dictionary mapping leaf node offsets to tuples of
                `(data_offset,data_block_size)`
 
        
        In addition, if `memorize` is set to `True,` the following dictionaries
        are additionally populated to speed subsequent searches:
 
            `self.node_genome_boundaries`
                Dictionary all node offsets to tuples of
                `(start_chromosome_id,start_genomic_base,end_chromosome_id,end_genomic_base)`
    
            
            `self.node_child_offsets`
                Dictionary mapping non-leaf node offsets to lists of offsets 
                of their child nodes
            
            `self.node_parent_offsets`
                Dictionary mapping all node offsets to the offset of the parent node
                for reverse tree traversal
        
        
        Parameters
        ----------
        start_node_offset : int, optional
            Start offset of first node block (Default: top of |RTree|) 
        """
        if start_node_offset is None:
            start_node_offset = self.tree_offset
        
        # otherwise, search the tree
        fh = open(self.filename,"rb")
        fh.seek(start_node_offset)
        node_info = RTreeNodeFormatFactory(fh,self._byte_order)
        
        for _ in range(node_info["count"]):
            current_offset = fh.tell()
            
            if node_info["is_leaf"] == True:
                # save leaf info
                new_node = RTreeLeafFactory(fh)
                self.leaf_data[current_offset] = (new_node["data_offset"],
                                                  new_node["data_size"])
            else:
                # if saving data, forward-link to child
                new_node = RTreeNonLeafFactory(fh)
                child_offset = new_node["child_data_offset"]
                if self.memorize == True:
                    try:
                        self.node_child_offsets[current_offset].append(child_offset)
                    except KeyError:
                        self.node_child_offsets[current_offset] = [child_offset]
                
                # recurse on child
                self.find_leaves(start_node_offset=child_offset)

            if self.memorize == True:
                # backlink to parent
                self.node_parent_offsets[current_offset] = start_node_offset
                
                # save genome boundaries
                self.node_genome_boundaries[current_offset] = (new_node["start_chrom_id"],
                                                               new_node["start_base"],
                                                               new_node["end_chrom_id"],
                                                               new_node["end_base"]
                                                               )                    
        
        fh.close()

#===============================================================================
# INDEX: Factories for various record formats
#===============================================================================

HeaderFactory = BinaryParserFactory("BigBedHeader","IHH3QHHQQIQ",["magic",
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
                                                                  "reserved"
                                                                  ])
"""Parse headers for BigWig and `BigBed`_ files"""

ZoomHeaderFactory = BinaryParserFactory("ZoomHeader","2H2Q",["reduction_level",
                                                             "reserved",
                                                             "data_offset",
                                                             "index_offset"])
"""Parse Zoomlevel headers in BigWig files"""


TotalSummaryFactory = BinaryParserFactory("TotalSummary","5Q",["bases_covered",
                                                               "min_val",
                                                               "max_val",
                                                               "sum_data",
                                                               "sum_squares"])
"""Parse 'Total Summary' tables in BigWig and `BigBed`_ files"""

BPlusTreeHeaderFactory = BinaryParserFactory("BPlusTreeHeader",
                                            "4I2Q",
                                            ["magic",
                                             "block_size",
                                             "key_size",
                                             "val_size",
                                             "num_chroms",
                                             "reserved"])
"""Parse headers of |BPlusTree|"""

BPlusTreeNodeFormatFactory  = BinaryParserFactory("BPlusTreeNodeFormat",
                                                  "?BH",
                                                  ["is_leaf",
                                                   "reserved",
                                                   "child_count"])
"""Determine format of BPlus Tree nodes in a datablock of a `BigBed`_ or BigWig file"""

RTreeHeaderFactory = BinaryParserFactory("RTreeHeader",
                                         "IIQ4IQII",
                                        ["magic",
                                         "block_size",
                                         "num_nodes",
                                         "start_chrom_id",
                                         "start_base",
                                         "end_chrom_id",
                                         "end_base",
                                         "end_file_offset",
                                         "items_per_leaf",
                                         "reserved"])
"""Parse headers of |RTree|"""

RTreeNodeFormatFactory = BinaryParserFactory("RTreeNodeFormat",
                                            "BBH",
                                            ["is_leaf","reserved","count"]
                                            )
"""Determine format of R Tree nodes in a data block of a `BigBed`_ or BigWig file"""


RTreeLeafFactory = BinaryParserFactory("RTreeLeaf",
                                       "4I2Q",
                                      ["start_chrom_id",
                                       "start_base",
                                       "end_chrom_id",
                                       "end_base",
                                       "data_offset",
                                       "data_size"
                                       ])
"""Parse leaves of |RTree|"""

RTreeNonLeafFactory = BinaryParserFactory("RTreeNonLeaf",
                                          "4IQ",
                                          ["start_chrom_id",
                                           "start_base",
                                           "end_chrom_id",
                                           "end_base",
                                           "child_data_offset"])
"""Parse non-leaf nodes of |RTree|"""


