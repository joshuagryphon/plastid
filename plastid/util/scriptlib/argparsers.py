#!/usr/bin/env python
"""This module contains pairs of functions that enable command-line scripts
to open and process various types of data. For each argument parser, one 
function creates the parser; the other parses it into an object.

===========================================================   ======================================  =======================================   ==========================================
**Capability**                                                **Parser building function**            **Argument parsing function**             **Return type**
-----------------------------------------------------------   --------------------------------------  ---------------------------------------   ------------------------------------------------------------ 
Open read alignments or count files                           :func:`get_alignment_file_parser`       :func:`get_genome_array_from_args`        |GenomeArray|, |BAMGenomeArray|, or |SparseGenomeArray|
(`BAM`_, `Bowtie`_, `bedGraph`_, `Wiggle`_ formats)
-----------------------------------------------------------   --------------------------------------  ---------------------------------------   ------------------------------------------------------------ 
Open and assemble transcript models                           :func:`get_annotation_file_parser`      :func:`get_transcripts_from_args`         iterator over |Transcripts|
(`BED`_, `BigBed`_, `GTF2`_, `GFF3`_, `PSL`_)
-----------------------------------------------------------   --------------------------------------  ---------------------------------------   ------------------------------------------------------------ 
Open and assemble arbitrary genomic features                  :func:`get_segmentchain_file_parser`    :func:`get_segmentchains_from_args`       iterator over |SegmentChains|
(`BED`_, `BigBed`_, `GTF2`_, `GFF3`_, `PSL`_)
-----------------------------------------------------------   --------------------------------------  ---------------------------------------   ------------------------------------------------------------ 
Open sequence files                                           :func:`get_sequence_file_parser`        :func:`get_seqdict_from_args`             dict-like of :class:`Bio.SeqRecord.SeqRecord` objects
(`fasta`_, `fastq`_, `twobit`_, genbank, embl)
-----------------------------------------------------------   --------------------------------------  ---------------------------------------   ------------------------------------------------------------ 
Hash genomic features by position                             :func:`get_mask_file_parser`            :func:`get_mask_hash_from_args`           |GenomeHash|, |TabixGenomeHash|, or |BigBedGenomeHash|
(`BED`_, `BigBed`_, `GTF2`_, `GFF3`_, `PSL`_)
===========================================================   ======================================  =======================================   ==========================================


 .. TODO remove this passage or decide on formatting

    Capabilities provided include:

        Importing read alignments or counts
            :func:`get_alignment_file_parser`
                Get :class:`~argparse.ArgumentParser` that opens alignment or
                count data in `BAM`_, `Bowtie`_, `Wiggle`_, and `bedGraph`_ files
        
            :func:`get_genome_array_from_args`
                Return a |GenomeArray|, |SparseGenomeArray|, or |BAMGenomeArray|,
                from command-line arguments parsed by a parser created with 
                :func:`get_alignment_file_parser`,

        Importing transcript models
            :func:`get_annotation_file_parser`
                Get :class:`~argparse.ArgumentParser` for opening transcript models
                from `GTF2`_, `GFF3`_, `BigBed`_, and `BED`_ format annotation files
                
            :func:`get_transcripts_from_args`
                Return |Transcripts| from command-line arguments parsed
                by a parser created with :func:`get_annotation_file_parser`

        Importing arbitrary regions of interest
            :func:`get_segmentchain_file_parser`
                Get :class:`~argparse.ArgumentParser` for opening regions of interest
                as |SegmentChains| from `GTF2`_, `GFF3`_, `BigBed`_, `BED`_, and `PSL`_
                format annotation files

            :func:`get_segmentchains_from_args`
                Return |SegmentChains| from arguments parsed by a parser created with
                :func:`get_segmentchain_file_parser`

        Opening sequence files
            :func:`get_sequence_file_parser`
                Get :class:`~argparse.ArgumentParser` for opening files of
                nucleic acid sequence

            :func:`get_seqdict_from_args`
                Return a dict-like object of SeqRecords of sequences from 
                a parser mad with :func:`get_sequence_file_parser`


        Hashing regions of the genome that should be excluded from analyses
            :func:`get_mask_file_parser`
                Get :class:`~argparse.ArgumentParser` to open a :term:`mask file` 
                in `GTF2`_, `GFF3`_, `BigBed`_, `BED`_, or `PSL`_ format.
                 
            :func:`get_genome_hash_from_mask_args`
                Return a |GenomeHash|, |TabixGenomeHash|, or |BigBedGenomeHash|
                of genomic regions to mask from analyses from a parser made with
                :func:`get_mask_file_parser`

Example
-------     
To add any these capabilities to your command line scripts, follow the steps
below. In this example, we create a script that can fetch |Transcripts| from
all of the annotation file formats supported by :data:`plastid`:

  #. Import one of the function pairs above::

         import argparse
         from plastid.util.scriptlib.argparsers import get_annotation_file_parser, get_transcripts_from_args


  #. Use the first function to create an  :class:`~argparse.ArgumentParser`,
     and supply this object as a `parent` when you build your script's
     :py:class:`~argparse.ArgumentParser`::

         # create annotation file parser
         annotation_file_parser = get_annotation_file_parser(disabled=["some_option"])
        
         # create my own parser, incorporating flags from annotation_file_parser
         my_own_parser = argparse.ArgumentParser(parents=[annotation_file_parser])     
    
         # add script-specific arguments
         my_own_parser.add_argument("--foo",type=int,default=5,help="Some option")
         my_own_parser.add_argument("--bar",type=str,default="a string",help="Another option")

     
  #. Then, use the second function to parse the arguments::

         # parse args
         args = parser.parse_args()
        
         # get transcript objects from arguments
         transcripts = get_transcripts_from_args(args)
        
         pass # rest of your script


Your script will then be able process various sorts  of count/alignment and annotation files
from the arguments.



See Also
--------
:py:mod:`argparse`
    Python documentation on argument parsing

:py:mod:`~plastid.genomics.genome_array`
    Data structures for import of read alignments or counts

:py:mod:`~plastid.genomics.roitools`
    Data structures describing transcripts or genomic regions of interest

:py:mod:`~plastid.genomics.genome_hash`
    Fetch feature annotations overlapping genomic regions of interest

:py:obj:`plastid.bin`
    Source code of command-line scripts, for further examples
"""
import argparse
import warnings
import sys
from collections import OrderedDict
from plastid.util.services.exceptions import MalformedFileError, ArgumentWarning
#from plastid.genomics.roitools import SegmentChain, Transcript
from plastid.util.io.openers import opener, NullWriter


#===============================================================================
# INDEX: String constants used in parsers below
#===============================================================================

_MAPPING_RULE_DESCRIPTION = """For BAM or bowtie files, one of the mutually exclusive read mapping choices
(`fiveprime_variable`, `fiveprime`, `threeprime`, or `center`) is required.
`--offset`, `--nibble`, `--min_length` and `--max_length` are optional."""

_DEFAULT_ALIGNMENT_FILE_PARSER_DESCRIPTION = "Open alignment or count files and optionally set mapping rules"
_DEFAULT_ALIGNMENT_FILE_PARSER_TITLE = "alignment mapping rules (for BAM & bowtie files)"

_DEFAULT_ANNOTATION_PARSER_DESCRIPTION = "Open one or more genome annotation files"
_DEFAULT_ANNOTATION_PARSER_TITLE = "annotation file options (one or more annotation files required)"

_MASK_PARSER_TITLE = "mask file options (optional)"
_MASK_PARSER_DESCRIPTION = """Add mask file(s) that annotate regions that should be excluded from analyses
(e.g. repetitive genomic regions)."""

_DEFAULT_SEQUENCE_PARSER_TITLE = "sequence options"
_DEFAULT_SEQUENCE_PARSER_DESCRIPTION = ""


GFF_SORT_MESSAGE = """Sort and index your GTF2/GFF with Tabix as follows:

    $ sort -k1,1 -k4,4n my_file.FORMAT | bgzip > my_file_sorted.FORMAT.gz
    $ tabix -p gff my_file_sorted.FORMAT.gz

See http://www.htslib.org/doc/tabix.html for download and documentation of tabix and bgzip."""


#===============================================================================
# INDEX: Alignment/count file parser, and helper functions
#===============================================================================


def get_alignment_file_parser(input_choices=("BAM","bowtie","wiggle"),
                              disabled=[],
                              prefix="",
                              title=_DEFAULT_ALIGNMENT_FILE_PARSER_TITLE,
                              description=_DEFAULT_ALIGNMENT_FILE_PARSER_DESCRIPTION,
                              return_subparsers=False):
    """Return an :py:class:`~argparse.ArgumentParser` that opens
    alignment (`BAM`_ or `bowtie`_) or count (`Wiggle`_, `bedGraph`_) files.
     
    In the case of `bowtie`_ or `BAM`_ import, parse arguments for mapping rules
    (e.g. fiveprime end mapping, threeprime end mapping, et c) and optional 
    read length filters
    
    Parameters
    ----------
    input_choices : list, optional
        list of permitted alignment file type choices
        
    disabled : list, optional
        list of parameter names that should be disabled from parser,
        without preceding dashes

    prefix : str, optional
        string prefix to add to default argument options (Default: "")
    
    title : str, optional
        title for option group (used in command-line help screen)
            
    description : str, optional
        description of parser (used in command-line help screen)
        
    return_subparsers : bool, optional
        if True, additionally return a dictionary of subparser option groups,
        to which additional options may be added (Default: False)
            
    Returns
    -------
    :class:`argparse.ArgumentParser`
    
    
    See also
    --------
    get_genome_array_from_args
        function that parses the :py:class:`~argparse.Namespace`
        returned by this :py:class:`~argparse.ArgumentParser`
    """
    alignment_file_parser = argparse.ArgumentParser(description=description,
                                                    add_help=False)

    subparsers = { }
    if len({ "BAM", "bowtie" } & set(input_choices)) > 0:
        subparsers["mapping"] = alignment_file_parser.add_argument_group(title=title,
                                                           description=_MAPPING_RULE_DESCRIPTION)
        
        # dictionary of options for mapping. Defined initially as a dict
        # so we can disable anything programmatically above
        map_option_dict = OrderedDict([
                            ("fiveprime_variable" , dict(action="store_const",
                                                        const="fiveprime_variable",
                                                        dest="%smapping" % prefix,
                                                        help="Map read alignment to a variable offset from 5' position of read, "+
                                                             "with offset determined by read length. Requires `--offset` below")),
                            ("fiveprime"        , dict(action="store_const",
                                                        const="fiveprime",
                                                        dest="%smapping" % prefix,
                                                        help="Map read alignment to 5' position.")),
                            ("threeprime"       , dict(action="store_const",
                                                        const="threeprime",
                                                        dest="%smapping" % prefix,
                                                        help="Map read alignment to 3' position")),
                            ("center"           , dict(action="store_const",
                                                        const="center",
                                                        dest="%smapping" % prefix,
                                                        help="Subtract N positions from each end of read, "+
                                                             "and add 1/(length-N), to each remaining position, "+
                                                             "where N is specified by `--nibble`")),
                            ("offset"           , dict(default=0,
                                                        metavar="OFFSET",
                                                        help="For `--fiveprime` or `--threeprime`, provide an integer "+
                                                          "representing the offset into the read, "+
                                                          "starting from either the 5\' or 3\' end, at which data "+
                                                          "should be plotted. For `--fiveprime_variable`, "+
                                                          "provide the filename of a two-column tab-delimited text "+
                                                          "file, in which first column represents read length or the "+
                                                          "special keyword `'default'`, and the second column represents "+
                                                          "the offset from the five prime end of that read length at which the read should be mapped.")),
                            ("nibble"           , dict(type=int,
                                                        default=0,
                                                        metavar="N",
                                                        help="For use with `--center` only. nt to remove from each "+
                                                             "end of read before mapping (Default: 0)")),
                            ("min_length"       , dict(type=int,
                                                       default=25,
                                                       metavar="N",
                                                       help="Minimum read length required to be included"+
                                                            " (Default: 25)")),
                            ("max_length"       , dict(type=int,
                                                       default=100,
                                                       metavar="N",
                                                       help="Maximum read length permitted to be included"+
                                                            " (Default: 100)")),
                            ])

        for k,v in filter(lambda x: x[0] not in disabled,map_option_dict.items()):
            subparsers["mapping"].add_argument("--%s%s" % (prefix,k),**v)

    alignment_option_dict = OrderedDict([
                        ("count_files"     , dict(type=str,
                                                  default=[],
                                                  nargs="+",
                                                  help="One or more count or alignment file(s) from a single sample or set of samples to be pooled.")),
                        ("countfile_format", dict(choices=input_choices,
                                                  default="BAM",
                                                  help="Format of file containing alignments or counts (default: BAM)")),
                        ("big_genome"       , dict(action="store_true",
                                                   default=False,
                                                   help="Use slower but memory-efficient implementation "+
                                                        "for big genomes (e.g. >20 megabases; irrelevant "+
                                                        "for BAM files), or for memory-limited computers")),
                        ("normalize"        , dict(action="store_true",
                                                   help="Whether counts should be normalized"+
                                                        " to counts per million (usually not. default: False)",
                                                   default=False)),
                        ])

    for k,v in filter(lambda x: x[0] not in disabled,alignment_option_dict.items()):
        alignment_file_parser.add_argument("--%s%s" % (prefix,k),**v)

    if return_subparsers == True:
        return alignment_file_parser, subparsers
    else:
        return alignment_file_parser


def get_genome_array_from_args(args,prefix="",disabled=[],printer=NullWriter()):
    """Return a |GenomeArray|, |SparseGenomeArray| or |BAMGenomeArray|
    from arguments parsed by :py:func:`get_alignment_file_parser`
    
    Parameters
    ----------
    args : :py:class:`argparse.Namespace`
        Namespace object from :py:func:`get_alignment_file_parser`

    prefix : str, optional
        string prefix to add to default argument options (Default: "")
        Must be same prefix that was added in call to :py:func:`get_alignment_file_parser`
        (Default: "")

    disabled : list, optional
        list of parameter names that were disabled when the argparser was created
        in :py:func:`get_alignment_file_parser`. (Default: ``[]``)
        
    printer : file-like, optional
        A stream to which stderr-like info can be written (default: |NullWriter|) 
    
    
    Returns
    -------
    |GenomeArray|, |SparseGenomeArray|, or |BAMGenomeArray|
    
    
    See Also
    --------
    get_alignment_file_parser
        Function that creates :py:class:`~argparse.ArgumentParser` whose output
        :py:class:`~argparse.Namespace` is processed by this function        
    """
    args = PrefixNamespaceWrapper(args,prefix)
    import pysam
    from plastid.util.io.filters import CommentReader
    from plastid.genomics.genome_array import GenomeArray, SparseGenomeArray,\
                                               BAMGenomeArray,\
                                               SizeFilterFactory, CenterMapFactory,\
                                               FivePrimeMapFactory, ThreePrimeMapFactory,\
                                               VariableFivePrimeMapFactory,\
                                               five_prime_map,  \
                                               three_prime_map, \
                                               center_map,      \
                                               variable_five_prime_map
    
    # require at least one countfile
    if len(args.count_files) == 0:
        printer.write("Please include at least one input file.")
        sys.exit(1)
    
    # require mapping rules unless wiggle
    if args.mapping is None and args.countfile_format != "wiggle":
        printer.write("Please specify a read mapping rule.")
        sys.exit(1)
    
    if "countfile_format" not in disabled and args.countfile_format == "BAM":
        count_files = [pysam.Samfile(X,"rb") for X in args.count_files]
        try:
            gnd = BAMGenomeArray(count_files)
        except ValueError:
            printer.write("Input BAM file(s) not indexed. Please index via:")
            printer.write("")
            for fn in args.count_files:
                printer.write("    samtools index %s" % fn)
            printer.write("")
            printer.write("Exiting.")
            sys.exit(1)
            
        size_filter = SizeFilterFactory(min=args.min_length,max=args.max_length)
        gnd.add_filter("size:%s-%s" % (args.min_length,args.max_length) ,size_filter)
        if args.mapping == "fiveprime":
            map_function = FivePrimeMapFactory(int(args.offset))
        elif args.mapping == "threeprime":
            map_function = ThreePrimeMapFactory(int(args.offset))
        elif args.mapping == "center":
            map_function = CenterMapFactory(args.nibble)
        elif args.mapping == "fiveprime_variable":
            if str(args.offset) == "0":
                printer.write("Please specify a filename to use for fiveprime variable offsets in --offset.")
                sys.exit(1)
            offset_dict = _parse_variable_offset_file(CommentReader(open(args.offset)))
            map_function = VariableFivePrimeMapFactory(offset_dict)
        gnd.set_mapping(map_function)
        
        
    else:
        if "big_genome" not in disabled and args.big_genome == True:
            gnd = SparseGenomeArray()
        else:
            gnd = GenomeArray()
            
        if "countfile_format" not in disabled and args.countfile_format == "wiggle":
            for align_file in args.count_files:
                printer.write("Opening wiggle files %s..." % align_file)
                with open("%s_fw.wig" % align_file) as fh:
                    gnd.add_from_wiggle(fh,"+")
                with open("%s_rc.wig" % align_file) as fh:
                    gnd.add_from_wiggle(fh,"-")
        else:
            trans_args = { "nibble" : int(args.nibble) }
            if args.mapping == "fiveprime_variable":
                transformation = variable_five_prime_map
                if str(args.offset) == "0":
                    printer.write("Please specify a filename to use for fiveprime variable offsets in --offset.")
                    sys.exit(1)
                else:
                    with open(args.offset) as myfile:
                        trans_args["offset"] = _parse_variable_offset_file(CommentReader(myfile))
            else:
                trans_args["offset"] = int(args.offset)
                if args.mapping == "fiveprime":
                    transformation = five_prime_map
                elif args.mapping == "threeprime":
                    transformation = three_prime_map
                elif args.mapping == "entire":
                    transformation = center_map
                elif args.mapping == "center":
                    transformation = center_map
                else:
                    transformation = center_map
        
            for infile in args.count_files:
                with opener(infile) as my_file:
                    if args.countfile_format == "bowtie":
                        gnd.add_from_bowtie(my_file,transformation,min_length=args.min_length,max_length=args.max_length,**trans_args)
        
        printer.write("Counted %s total reads..." % gnd.sum())
        
    if "normalize" not in disabled and args.normalize == True:
        printer.write("Normalizing to reads per million...")
        gnd.set_normalize(True)
    
    return gnd


#===============================================================================
# INDEX: Annotation file parser, and helper functions
#===============================================================================


def get_annotation_file_parser(input_choices=["BED","BigBed","GTF2","GFF3"],
                               disabled=[],
                               prefix="",
                               title=_DEFAULT_ANNOTATION_PARSER_TITLE,
                               description=_DEFAULT_ANNOTATION_PARSER_DESCRIPTION,
                               return_subparsers=False):
    """Return an :py:class:`~argparse.ArgumentParser` that opens
    annotation files from `BED`_, `BigBed`_, `GTF2`_, or `GFF3`_ formats
     
    Parameters
    ----------
    input_choices : list, optional
        list of permitted alignment file type choices.
        (Default: '["BED","BigBed","GTF2","GFF3"]'). 'PSL'_ may also be added
        
    disabled : list, optional
        list of parameter names that should be disabled from parser
        without preceding dashes

    prefix : str, optional
        string prefix to add to default argument options (Default: `''`)
    
    title : str, optional
        title for option group (used in command-line help screen)
        
    description : str, optional
        description of parser (used in command-line help screen)
    
    return_subparsers : bool, optional
        if True, additionally return a dictionary of subparser option groups,
        to which additional options may be added (Default: `False`)
    
    Returns
    -------
    :class:`argparse.ArgumentParser`
    
    
    See also
    --------
    get_transcripts_from_args
        function that parses the :py:class:`~argparse.Namespace` returned
        by this :py:class:`~argparse.ArgumentParser`
    """
    from plastid.readers.gff import _DEFAULT_GFF3_GENE_TYPES,\
                                 _DEFAULT_GFF3_TRANSCRIPT_TYPES,\
                                 _DEFAULT_GFF3_EXON_TYPES,\
                                 _DEFAULT_GFF3_CDS_TYPES

    annotation_file_parser = argparse.ArgumentParser(add_help=False)
    """Open genome annotation files in %s format""" % ", ".join(input_choices)
    
    subparsers = { "annotation" : annotation_file_parser.add_argument_group(title=title,description=description) }

    option_dict = OrderedDict([
                    ("annotation_files"     , dict(metavar="infile.[%s]" % " | ".join(input_choices),# | psl]",
                                                  type=str,nargs="+",default=[],
                                                  help="Zero or more annotation files (max 1 file if BigBed)")),                               
                    ("annotation_format"   , dict(choices=input_choices,
                                                  default="GTF2",
                                                  help="Format of %sannotation_files (default: GTF2). Note: GFF3 assembly assumes SO v.2.5.2 feature ontologies, which may or may not match your specific file." % prefix)),    
                    ("add_three"           , dict(default=False,
                                                  action="store_true",
                                                  help="If supplied, coding regions will be extended by 3 nucleotides at their 3\' ends (except for GTF2 files that explicitly include `stop_codon` features). Use if your annotation file excludes stop codons from CDS.")),
                    ("tabix"               , dict(default=False,
                                                  action="store_true",
                                                  help="%sannotation_files are tabix-compressed and indexed (Default: False). Ignored for BigBed files." % prefix)),
                    ("sorted"              , dict(default=False,
                                                  action="store_true",
                                                  help="%sannotation_files are sorted by chromosomal position (Default: False)" % prefix))
                   ])

    if "BED" in input_choices:
        option_dict["bed_extra_columns"] = dict(default=0,
                                                nargs="+",
                                                help="Number of extra columns in BED file (e.g. in custom ENCODE formats) "+
                                                     "or list of names for those columns. (Default: 0).")
    if "GFF3" in input_choices:
        option_dict["gff_transcript_types"] = dict(type=str,
                                                  default=_DEFAULT_GFF3_TRANSCRIPT_TYPES,
                                                  nargs="+",
                                                  help="GFF3 feature types to include as transcripts, even "+\
                                                 "if no exons are present (for GFF3 only; default: use SO v2.5.3 specification)")
        option_dict["gff_exon_types"] =       dict(type=str,
                                                  default=_DEFAULT_GFF3_EXON_TYPES,
                                                  nargs="+",
                                                  help="GFF3 feature types to include as exons (for GFF3 only; default: use SO v2.5.3 specification)")
        option_dict["gff_cds_types"] =       dict(type=str,
                                                  default=_DEFAULT_GFF3_CDS_TYPES,
                                                  nargs="+",
                                                  help="GFF3 feature types to include as CDS (for GFF3 only; default: use SO v2.5.3 specification)")
        
    for k,v in filter(lambda x: x[0] not in disabled,option_dict.items()):
        subparsers["annotation"].add_argument("--%s%s" % (prefix,k),**v)
    
    if return_subparsers == True:
        return annotation_file_parser, subparsers
    else:
        return annotation_file_parser

def get_transcripts_from_args(args,prefix="",disabled=[],printer=NullWriter(),return_type=None,require_sort=False):
    """Return a list of |Transcript| objects from arguments parsed by :py:func:`get_annotation_file_parser`
    
    Parameters
    ----------
    args : :py:class:`argparse.Namespace`
        Namespace object from :py:func:`get_annotation_file_parser`
    
    prefix : str, optional
        string prefix to add to default argument options.
        Must be same prefix that was added in call to :py:func:`get_annotation_file_parser`
        (Default: `''`)
        
    disabled : list, optional
        list of parameter names that were disabled when the annotation file
        parser was created by :py:func:`get_annotation_file_parser`. 
        (Default: `[]`)
            
    printer : file-like, optional
        A stream to which stderr-like info can be written (Default: |NullWriter|) 
    
    return_type : |SegmentChain| or subclass, optional
        Type of object to return (Default: |Transcript|)

    require_sort : bool, optional
        If True, quit if the annotation file(s) are not sorted or indexed
    
    Returns
    -------
    iterator
        |Transcript| objects, either in order of appearance (if input was a
        `BED`_, `BigBed`_, or `PSL`_ file), or sorted lexically by chromosome,
        start coordinate, end coordinate, and then strand (if input was `GTF2`_
        or `GFF3`_).
    
    
    See Also
    --------
    get_annotation_file_parser
        Function that creates :py:class:`argparse.ArgumentParser` whose output
        :py:class:`~argparse.Namespace` is processed by this function    
    """
    if return_type is None:
        from plastid.genomics.roitools import Transcript
        return_type = Transcript

    if prefix != "":
        args = PrefixNamespaceWrapper(args,prefix)

    if require_sort == True and 'sorted' not in disabled:
        if args.annotation_format in ("BED","GTF2","GFF3") and \
            args.sorted == False and 'tabix' not in disabled and\
            args.tabix == False:
            printer.write("Using unsorted/unindexed annotation files requires impractical amounts of memory.")
            if args.annotation_format == "BED":
                printer.write("""Convert BED to BigBed using Jim Kent's bedToBigBed utility as follows:

    $ sort -k1,1 -k2,2n my_file > my_file_sorted.bed
    $ bedToBigBed my_file_sorted.bed chrom.sizes my_file_sorted.bb

See https://github.com/ENCODE-DCC/kentUtils/tree/master/src/product/scripts
for download & documentation of Kent utilities""")
                sys.exit(1)
            else:
                printer.write(GFF_SORT_MESSAGE.replace("FORMAT",args.annotation_format.lower()))
                sys.exit(1)

    printer.write("Parsing features in %s..." % ", ".join(args.annotation_files))
    
    if "tabix" not in disabled:
        tabix = args.tabix
    else:
        tabix = False
        
    if "add_three" not in disabled:
        add_three = args.add_three
    else:
        add_three = False

    if "bed_extra_columns" not in disabled:
        bed_extra_columns = args.bed_extra_columns
        if not(isinstance(bed_extra_columns,list)):
            try:
                bed_extra_columns = int(bed_extra_columns)
            except ValueError:
                pass
    else:
        bed_extra_columns = 0
    
    if args.annotation_format.lower() == "bigbed":
        if len(args.annotation_files) > 1:
            printer.write("Bad arguments: we can only process one BigBed file.")
            sys.exit(2)
        if tabix == True:
            warnings.warn("Tabix compression is incompatible with BigBed files. Ignoring.",ArgumentWarning)

        from plastid.readers.bigbed import BigBedReader
        transcripts = BigBedReader(args.annotation_files[0],
                                   return_type=return_type,
                                   cache_depth=1,
                                   add_three_for_stop=add_three,
                                   printer=printer)
        
    elif tabix == True:
        import pysam
        #streams = [pysam.tabix_iterator(opener(X), lambda x,y: x) for X in args.annotation_files] # used to work in earlier pysam
        # string parsing by supplying None instead of `asTuple()` no longer works
        # nor do anonymous lambda functions
        streams = [pysam.tabix_iterator(opener(X), pysam.asTuple()) for X in args.annotation_files]
    else:
        streams = (opener(X) for X in args.annotation_files)

    if args.annotation_format in ("GFF3","GTF2"):
        from plastid.readers.gff import GFF3_TranscriptAssembler, GTF2_TranscriptAssembler
        if 'sorted' not in disabled and args.sorted == False and \
           'tabix' not in disabled and args.tabix == False:
            msg = """Transcript assembly on FORMAT files can require a lot of memory.
Consider using a sorted file with the '--sorted' flag and/or tabix-compression.
"""
            msg += GFF_SORT_MESSAGE.replace("FORMAT",args.annotation_format.lower())
            warnings.warn(msg,ArgumentWarning)
    
    if args.annotation_format.lower() == "gff3":
        transcripts = GFF3_TranscriptAssembler(*streams,
                                               transcript_types=args.gff_transcript_types,
                                               exon_types=args.gff_exon_types,
                                               cds_types=args.gff_cds_types,
                                               printer=printer,
                                               add_three_for_stop=add_three,
                                               tabix=tabix,
                                               return_type=return_type,
                                               is_sorted=args.sorted)
    elif args.annotation_format.lower() == "gtf2":
        transcripts = GTF2_TranscriptAssembler(*streams,
                                               printer=printer,
                                               tabix=tabix,
                                               return_type=return_type,
                                               add_three_for_stop=add_three,
                                               is_sorted=args.sorted)
        
    elif args.annotation_format.lower() == "bed":
        from plastid.readers.bed import BED_Reader
        transcripts = BED_Reader(*streams,
                                 add_three_for_stop=add_three,
                                 tabix=tabix,
                                 return_type=return_type,printer=printer,
                                 extra_columns=bed_extra_columns)

    elif args.annotation_format.lower() == "psl":
        from plastid.readers.psl import PSL_Reader
        transcripts = PSL_Reader(*streams,
                                 tabix=tabix,
                                 return_type=return_type,printer=printer)
        
    return transcripts

def get_segmentchain_file_parser(input_choices=["BED","BigBed","GTF2","GFF3","PSL"],
                                 disabled=[],
                                 prefix="",
                                 title=_DEFAULT_ANNOTATION_PARSER_TITLE,
                                 description=_DEFAULT_ANNOTATION_PARSER_DESCRIPTION):
    """Create an :class:`~argparse.ArgumentParser` to open annotation files as |SegmentChains|
    
    Parameters
    ----------
    input_choices : list, optional
        list of permitted alignment file type choices
        (Default: `["BED","BigBed","GTF2","GFF3", "PSL"]`)
        
    disabled : list, optional
        list of parameter names that should be disabled from parser
        without preceding dashes

    prefix : str, optional
        string prefix to add to default argument options (Default: `''`)
    
    title : str, optional
        title for option group (used in command-line help screen)
        
    description : str, optional
        description of parser (used in command-line help screen)
 

    Returns
    -------
    :class:`argparse.ArgumentParser`
    
    
    See Also
    --------
    get_segmentchains_from_args
        function that parses the :py:class:`~argparse.Namespace` returned
        by this :py:class:`~argparse.ArgumentParser`
    """
    disabled.append([prefix+"add_three"])
    return get_annotation_file_parser(input_choices=input_choices,
                                      prefix=prefix,
                                      title=title,
                                      disabled=disabled,
                                      description=description)

def get_segmentchains_from_args(args,prefix="",disabled=[],printer=NullWriter(),require_sort=False):
    """Return a list of |SegmentChain| objects from arguments parsed by an
    :py:class:`~argparse.ArgumentParser` created by :py:func:`get_segmentchain_file_parser`
    
    Parameters
    ----------
    args : :py:class:`argparse.Namespace`
        Namespace object from :py:func:`get_segmentchain_file_parser`

    prefix : str, optional
        string prefix to add to default argument options.
        Must be same prefix that was added in call to :py:func:`get_segmentchain_file_parser`
        (Default: "")
        
    disabled : list, optional
        list of parameter names that were disabled when the annotation file
        parser was created by :py:func:`get_segmentchain_file_parser`. 
        (Default: ``[]``)
                
    printer : file-like
        A stream to which stderr-like info can be written (default: |NullWriter|) 
    
    require_sort : bool, optional
        If True, quit if the annotation file(s) are not sorted or indexed

    
    Returns
    -------
    iterator
        sequence of |SegmentChain| objects, either in order of appearance
        (if input was a BED or PSL file), or sorted lexically by chromosome,
        start coordinate, end coordinate, and then strand (if input was) GTF or GFF
    
    
    See Also
    --------
    get_segmentchain_file_parser
        Function that creates :py:class:`argparse.ArgumentParser` whose output
        :py:class:`~argparse.Namespace` is processed by this function    
    """
    from plastid.genomics.roitools import SegmentChain
    disabled.append([prefix+"add_three"])
    return get_transcripts_from_args(args,
                                     prefix=prefix,
                                     disabled=disabled,
                                     printer=printer,
                                     return_type=SegmentChain,
                                     require_sort=require_sort)

def get_mask_file_parser(prefix="mask_",disabled=[]):
    """Create an :class:`~argparse.ArgumentParser` to open annotation files that describe regions of the genome to mask from analyses
    
    Parameters
    ----------
    prefix : str, optional
        Prefix to add to default argument options (Default: `'mask_'`)
        
    disabled : list, optional
        list of parameter names to disable from the mask file parser 
        (Default: `[]`. `add_three` is always disabled.)

    Returns
    -------
    argparse.ArgumentParser
    
    See Also
    --------
    get_genome_hash_from_mask_args
        function that parses the :py:class:`~argparse.Namespace` returned
        by this :py:class:`~argparse.ArgumentParser`    
    """
    disabled.append(prefix+"add_three")
    return get_segmentchain_file_parser(prefix=prefix,
                                        disabled=disabled,
                                        input_choices=["BED","GTF2","GFF3","BigBed","PSL"],
                                        title=_MASK_PARSER_TITLE,
                                        description=_MASK_PARSER_DESCRIPTION)


def get_genome_hash_from_mask_args(args,prefix="mask_",printer=NullWriter()):
    """Return a |GenomeHash| of regions from command-line arguments

    Parameters
    ----------
    args : :py:class:`argparse.Namespace`
        Namespace object from :py:func:`get_mask_file_parser`
    
    prefix : str, optional
        string prefix to add to default argument options.
        Must be same prefix that was added in call to :py:func:`get_mask_file_parser`
        (Default: "mask_")
    
    printer : file-like
        A stream to which stderr-like info can be written (default: |NullWriter|) 


    Returns
    -------
    |GenomeHash|
        Hashed data structure of masked genomic regions
        
    
    See Also
    --------
    get_mask_file_parser
        Function that creates :py:class:`argparse.ArgumentParser` whose output
        :py:class:`~argparse.Namespace` is processed by this function  
    """
    from plastid.genomics.genome_hash import GenomeHash, BigBedGenomeHash, TabixGenomeHash
    tmp = PrefixNamespaceWrapper(args,prefix)
    if len(tmp.annotation_files) > 0:
        printer.write("Opening mask annotation file(s) %s..." % ", ".join(tmp.annotation_files))
        if tmp.annotation_format in ("BED","GTF2","GFF3") and tmp.tabix == False:
            msg = """Unindexed mask files can require lots of memory in large (e.g. mammalian) genomes.
Consider converting to BigBed or using tabix to index your mask file."""
            warnings.warn(msg,ArgumentWarning)

        if len(tmp.annotation_files) > 0:
            if tmp.annotation_format.lower() == "bigbed":
                if len(tmp.annotation_files) > 1:
                    printer.write("Bad arguments: we can only process one BigBed file.")
                    sys.exit(2)
                return BigBedGenomeHash(tmp.annotation_files[0])
            elif tmp.tabix == True:
                return TabixGenomeHash(tmp.annotation_files,tmp.annotation_format,printer=printer)
            else:
                hash_ivcs = get_segmentchains_from_args(args,prefix=prefix,printer=printer)
                return GenomeHash(hash_ivcs)
    else:
        return GenomeHash()
    

#===============================================================================
# INDEX: Sequence file parser
#===============================================================================

def get_sequence_file_parser(input_choices=["fasta","fastq","twobit","genbank","embl"],
                               disabled=[],
                               prefix="",
                               title=_DEFAULT_SEQUENCE_PARSER_TITLE,
                               description=_DEFAULT_SEQUENCE_PARSER_DESCRIPTION):
    """Return an :py:class:`~argparse.ArgumentParser` that opens
    annotation files from `BED`_, `BigBed`_, `GTF2`_, or `GFF3`_ formats
     
    Parameters
    ----------
    input_choices : list, optional
        list of permitted sequence file type choices.
        (Default: '["FASTA","twobit","genbank","embl"]').
        
    disabled : list, optional
        list of parameter names that should be disabled from parser
        without preceding dashes

    prefix : str, optional
        string prefix to add to default argument options (Default: `''`)
    
    title : str, optional
        title for option group (used in command-line help screen)
        
    description : str, optional
        description of parser (used in command-line help screen)
    
   
    Returns
    -------
    :class:`argparse.ArgumentParser`
    
    
    See also
    --------
    get_seqdict_from_args
        function that parses the :py:class:`~argparse.Namespace` returned
        by this :py:class:`~argparse.ArgumentParser`
    """
    sequence_file_parser = argparse.ArgumentParser(add_help=False)
    """Open sequence file in %s formats""" % ", ".join(input_choices)

    option_dict = OrderedDict([
                    ("sequence_file"     , dict(metavar="infile.[%s]" % " | ".join(input_choices),
                                                type=str,
                                                help="A file of DNA sequence")),                               
                    ("sequence_format"   , dict(choices=input_choices,
                                                default="fasta",
                                                help="Format of %ssequence_file (default: fasta)." % prefix)),    
                   ])

    for k,v in filter(lambda x: x[0] not in disabled,option_dict.items()):
        sequence_file_parser.add_argument("--%s%s" % (prefix,k),**v)
    
    return sequence_file_parser

def get_seqdict_from_args(args,index=True,prefix="",printer=NullWriter()):
    """Retrieve a dictionary-like object of sequences

    Parameters
    ----------
    args : :py:class:`argparse.Namespace`
        Namespace object from :py:func:`get_sequence_file_parser`
    
    prefix : str, optional
        string prefix to add to default argument options.
        Must be same prefix that was added in call to :py:func:`get_sequence_file_parser`
        (Default: "")
   
    index : bool, optional
        If sequence format is anything other than twobit, open with
        lazily-evaluating :func:`Bio.SeqIO.index` instead of
        :func:`Bio.SeqIO.to_dict` (Default: `True`)
        
    printer : file-like
        A stream to which stderr-like info can be written (default: |NullWriter|) 

    Returns
    -------
    dict-like
        Dictionary-like object mapping chromosome names to
        :class:`Bio.SeqRecord.SeqRecord`-like objects
    """
    args = PrefixNamespaceWrapper(args,prefix)
    printer.write("Opening sequence file '%s'." % args.sequence_file)
    if args.sequence_format == "twobit":
        from plastid.genomics.seqtools import TwoBitSeqRecordAdaptor
        return TwoBitSeqRecordAdaptor(args.sequence_file)
    else:
        from Bio import SeqIO
        if index == True:
            return SeqIO.index(args.sequence_file,args.sequence_format)
        else:
            return SeqIO.to_dict(SeqIO.parse(args.sequence_file,args.sequence_format))

#===============================================================================
# INDEX: Utility classes
#===============================================================================

class PrefixNamespaceWrapper(object):
    """Wrapper class to facilitate processing of :py:class:`~argparse.Namespace`
    objects created by :py:func:`get_alignment_file_parser` or
    :py:func:`get_annotation_file_parser` with non-empty ``prefix`` values,
    as if no prefix had been used.
    
    Attributes
    ----------
    namespace : :py:class:`~argparse.Namespace`
        Result of calling :py:meth:`argparse.ArgumentParser.parse_args`
    
    prefix : str
        Prefix that will be prepended to names of attributes of `self.namespace`
        before they are fetched. Must match prefix that was used in creation
        of the :py:class:`argparse.ArgumentParser` that created `self.namespace`
    
    See Also
    --------
    get_annotation_file_parser
    
    get_alignment_file_parser
    
    get_genome_array_from_args
    
    get_transcripts_from_args
    """
    
    def __init__(self,namespace,prefix):
        """Create a |PrefixNamespaceWrapper|
        
        Parameters
        ----------
        namespace : :py:class:`~argparse.Namespace`
            Result of calling :py:meth:`argparse.ArgumentParser.parse_args`
        
        prefix : str
            Prefix that will be prepended to items from the :py:class:`~argparse.Namespace`
            before they are checked 
        """
        self.namespace = namespace
        self.prefix = prefix
    
    def __getattr__(self,k):
        """Fetch an attribute from `self.namespace`, appending `self.prefix` to `k`
        before fetching
        
        Parameters
        ----------
        k : str
            Attribute to fetch
        """
        return getattr(self.namespace,"%s%s" % (self.prefix,k))


#===============================================================================
# INDEX: Utility functions
#===============================================================================


def _parse_variable_offset_file(fh):
    """Read a variable-offset text file into a dictionary.
    These text files contain two columns and are tab-delimited. The first column
    specifies the read length, or contains the special value `'default'`. The
    second column specifies the offset from the 5' end of that read length to 
    use.
    
    Parameters
    ----------
    fh : file-like
        open filehandle pointing to data
    
    Returns
    -------
    dict
        dictionary mapping sequencing read lengths to their 5' offsets
    """
    my_dict = {}
    for line in fh:
        if line.startswith("length"):
            continue
        items = line.strip("\n").split("\t")
        if len(items) != 2:
            name = getattr(fh,"__name__","Variable offset file")
            raise MalformedFileError(name,"More or fewer than two columns on line:\n\t%s" % line.strip("\n"))
        if items[0] == "length":
            continue
        key = items[0]
        try:
            key = key if key == "default" else int(key)
        except ValueError:
            name = getattr(fh,"__name__","Variable offset file")
            raise MalformedFileError(name,"Non integer value for key '%s' on line:\n\t%s" % (key,line.strip("\n")))
        if key in my_dict:
            name = getattr(fh,"__name__","Variable offset file")
            raise MalformedFileError(name,"multiple offsets defined for read length %s" % key)
        else:
            try:
                my_dict[key] = int(items[1])
            except ValueError:
                name = getattr(fh,"__name__","Variable offset file")
                raise MalformedFileError(name,"Non integer value for value '%s' on line:\n\t%s" % (items[1],line.strip("\n")))
            
    return my_dict
