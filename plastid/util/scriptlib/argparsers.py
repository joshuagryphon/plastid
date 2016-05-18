#!/usr/bin/env python
"""This module contains classes that:

  - build :class:`argparse.ArgumentParser` objects for various data types
    used in genomics
  
  - parse those arguments into useful file types
  

Arguments are grouped into the following sets:

    ===========================================================   ======================================
    **Parameter/argument set**                                    **Parser building class**            
    -----------------------------------------------------------   --------------------------------------
    Generic parameters (e.g. for error reporting, logging)        :class:`BaseParser`
    
    :term:`Read alignments` or :term:`count files <count file>`   :class:`AlignmentParser`
    
    Genomic feature or mask annotations                           :class:`AnnotationParser`
    
    Genomic sequence files                                        :class:`SequenceParser`
    
    Plotting parameters for charts                                :class:`PlottingParser`
    ===========================================================   ======================================


Example
-------     
To use any of these in your own command line scripts, follow these steps:

  #. Import one or more of the classes above::

         >>> import argparse
         >>> from plastid.util.scriptlib.argparsers import AnnotationParser


  #. Use the first function to create an  :class:`~argparse.ArgumentParser`,
     and supply this object as a `parent` when you build your script's
     :py:class:`~argparse.ArgumentParser`::

         >>> ap = AnnotationParser()
        
         # create annotation file parser
         >>> annotation_file_parser = ap.get_parser(disabled=["some_option_to_disable"])
        
         # create my own parser, incorporating flags from annotation_file_parser
         >>>> my_own_parser = argparse.ArgumentParser(parents=[annotation_file_parser])     
    
         # add script-specific arguments
         >>> my_own_parer.add_argument("positional_argument",type=str)
         >>> my_own_parser.add_argument("--foo",type=int,default=5,help="Some option")
         >>> my_own_parser.add_argument("--bar",type=str,default="a string",help="Another option")
     
  #. Then, use the second parse the arguments::

         >>> args = parser.parse_args()
        
         # get transcript objects from arguments
         # this will be an iterator over |Transcripts|
         >>> transcripts = ap.get_transcripts_from_args(args)
        
         >>> pass # rest of your script


Your script will then be able process whatever sorts of annotation files that
plastid currently supports.


See Also
--------
:py:mod:`argparse`
    Python documentation on argument parsing

:py:obj:`plastid.bin`
    Source code of command-line scripts, for further examples
"""
import sys
import os
import functools
import argparse
import warnings
import pkg_resources
import pysam

from plastid.util.services.exceptions import MalformedFileError, ArgumentWarning,\
                                             DataWarning, FileFormatWarning, filterwarnings
from plastid.util.services.decorators import deprecated
#from plastid.genomics.roitools import SegmentChain, Transcript
from plastid.util.io.openers import opener, NullWriter
from plastid.util.io.filters import CommentReader
from plastid.genomics.roitools import Transcript, SegmentChain
from plastid.readers.gff import _DEFAULT_GFF3_TRANSCRIPT_TYPES,\
                                _DEFAULT_GFF3_EXON_TYPES,\
                                _DEFAULT_GFF3_CDS_TYPES
                                
#===============================================================================
# INDEX: Constants used in parsers below
#===============================================================================

# Size above which we recommend sorting a GFF/GTF2 file
_GFF_SORT_SIZE = 100 * 1024 * 1024

_MAPPING_RULE_TITLE = "alignment mapping functions (BAM & bowtie files only)"
_MAPPING_RULE_DESCRIPTION = \
"""For BAM or bowtie files, one of the mutually exclusive read mapping functions
is required:
"""

_MAPPING_OPTION_TITLE = "filtering and alignment mapping options"
_MAPPING_OPTION_DESCRIPTION = \
"""
The remaining arguments are optional and affect the behavior of specific
mapping functions:
"""


_DEFAULT_ALIGNMENT_FILE_PARSER_DESCRIPTION = \
"Open alignment or count files and optionally set mapping rules"

_DEFAULT_ALIGNMENT_FILE_PARSER_TITLE = "count & alignment file options"



_DEFAULT_ANNOTATION_PARSER_DESCRIPTION = \
"Open one or more genome annotation files"

_DEFAULT_ANNOTATION_PARSER_TITLE = \
"annotation file options (one or more annotation files required)"

GFF_SORT_MESSAGE = \
"""Sort and index your GTF2/GFF with Tabix as follows:

    $ sort -k1,1 -k4,4n my_file.FORMAT | bgzip > my_file_sorted.FORMAT.gz
    $ tabix -p gff my_file_sorted.FORMAT.gz

See http://www.htslib.org/doc/tabix.html for download and documentation of
tabix and bgzip."""



_MASK_PARSER_TITLE = "mask file options (optional)"
_MASK_PARSER_DESCRIPTION = \
"""Add mask file(s) that annotate regions that should be excluded from analyses
(e.g. repetitive genomic regions)."""



_DEFAULT_SEQUENCE_PARSER_TITLE = "sequence options"
_DEFAULT_SEQUENCE_PARSER_DESCRIPTION = ""



_DEFAULT_PLOTTING_TITLE = "Plotting options"




#===============================================================================
# INDEX: Base class for parsers
#===============================================================================

class Parser(object):
    """Base class for argument parser factories used below
        
    Parameters
    ----------
    groupname : str, optional
        Name of argument group. If not `None`, an argument group with
        the specified name will be created and added to the parser.
        If not, arguments will be in the main group. 
        
    prefix : str, optional
        string prefix to add to default argument options (Default: "")

    disabled : list, optional
        list of parameter names that should be disabled from parser,
        without preceding dashes    
    """
    
    def __init__(self,groupname=None,prefix="",disabled=None,**kwargs):
        """Create a parser
        
        Parameters
        ----------
        groupname : str, optional
            Name of argument group. If not `None`, an argument group with
            the specified name will be created and added to the parser.
            If not, arguments will be in the main group. 
            
        prefix : str, optional
            string prefix to add to default argument options (Default: "")

        disabled : list, optional
            list of parameter names that should be disabled from parser,
            without preceding dashes
        """
        self.prefix = prefix
        self.disabled = [] if disabled is None else disabled
        self.groupname = groupname
        
        # define in __init__ of subclass
        self.arguments = []
    
    def get_parser(self,parser=None,groupname=None,arglist=None,title=None,description=None,**kwargs):
        """Create an populate :class:`argparse.ArgumentParser` with arguments
        
        Parameters
        ----------
        parser : :class:`argparse.ArgumentParser` or None, optional
            If `None`, a new parser will be created, and arguments will be added
            to it. If not `None`, arguments will be added to `parser`.
            (Default: `None`)
            
        groupname : str or None, optional
            If not `None`, default to `self.groupname`. If either `groupname`
            or `self.groupname` is not `None`, an option group with this name
            will be added to `parser`, and arguments added to that groupname
            instead of the main argument group of `parser`. In this case, `title`
            and `description` will be applied to the option group instead of to `parser`.
            Default : `None`)
            
        arglist : list, optional
            If not `None`, arguments in this list will be added to `parser`.
            Otherwise, arguments will be taken from `self.arguments`.
            
            The list should be a list of tuples of ('argument_name',dict_of_options),
            where `argument_name` is a string, and `dict_of_options` a dictionary
            of keyword arguments to pass to :meth:`argparse.ArgumentParser.add_argument`.
            
        title : str, optional
            Optional title for parser
            
        description : str, optional
            Optional description for parser


        kwargs : keyword arguments
            Additional arguments passed during creation of :class:`argparse.ArgumentParser`
        
        
        Returns
        -------
        :class:`argparse.ArgumentParser`
        """
        if groupname is None:
            groupname = self.groupname

        if parser is None:
            if groupname is None:
                parser = argparse.ArgumentParser(description=description,add_help=False,**kwargs)
            else:
                parser = argparse.ArgumentParser(add_help=False,**kwargs)
        
        addto = parser
        if groupname is not None:
            addto = parser.add_argument_group(title=title,description=description)
            
        arglist = self.arguments if arglist is None else arglist
        for arg_name, arg_opts in filter(lambda x: x[0] not in self.disabled,arglist):
            addto.add_argument("--%s%s" % (self.prefix,arg_name),**arg_opts)
            
        return parser

#===============================================================================
# INDEX: Alignment & count file parser
#===============================================================================


class AlignmentParser(Parser):
    """Parser for files containing read alignments or quantitative data.
    
    Checks for additional mapping rules and command-line arguments 
    by checking the entrypoints ``plastid.mapping_rules`` and
    ``plastid.mapping_options``
    
        
    Parameters
    ----------
    groupname : str, optional
        Name of argument group. If not `None`, an argument group with
        the specified name will be created and added to the parser.
        If not, arguments will be in the main group. 

    prefix : str, optional
        string prefix to add to default argument options (Default: "")

    disabled : list, optional
        list of parameter names that should be disabled from parser,
        without preceding dashes

    input_choices : list, optional
        list of permitted alignment file type choices for input
    
    allow_mapping : bool, optional
        Enable/disable user configuration of mapping rules (default: True)    
    """
    
    def __init__(self,prefix="",disabled=None,
                 input_choices=("BAM","bigwig","bowtie","wiggle"),
                 groupname="alignment_options",
                 allow_mapping=True):
        """Create a parser for read alignments and/or quantitative data
        
        Parameters
        ----------
        groupname : str, optional
            Name of argument group. If not `None`, an argument group with
            the specified name will be created and added to the parser.
            If not, arguments will be in the main group. 

        prefix : str, optional
            string prefix to add to default argument options (Default: "")

        disabled : list, optional
            list of parameter names that should be disabled from parser,
            without preceding dashes

        input_choices : list, optional
            list of permitted alignment file type choices for input
        
        allow_mapping : bool, optional
            Enable/disable user configuration of mapping rules (default: True)
        """
        Parser.__init__(self,groupname=groupname,prefix=prefix,disabled=disabled)
        self.input_choices = input_choices
        self.allow_mapping = allow_mapping
        
        self.bamfuncs = {}
        self.bowtiefuncs = {}
        
        self.arguments = [
            ("count_files"     , dict(type=str,
                                      default=[],
                                      nargs="+",
                                      help="One or more count or alignment file(s) from a single sample or set of samples to be pooled.")),
            ("countfile_format", dict(choices=input_choices,
                                      default="BAM",
                                      help="Format of file containing alignments or counts (Default: %(default)s)")),
            ("normalize"        , dict(action="store_true",
                                       help="Whether counts should be normalized"+
                                            " to counts per million (usually not. default: %(default)s)",
                                       default=False)),
            ("sum"              , dict(type=float,default=None,
                                       help="Sum used in normalization of counts and RPKM/RPNT calculations "+\
                                            "(Default: total mapped reads/counts in dataset)")),
                          
            ]
        
        length_ops = [
            ("min_length"       , dict(type=int,
                                       default=25,
                                       metavar="N",
                                       help="Minimum read length required to be included"+
                                            " (BAM & bowtie files only. Default: %(default)s)")),
            ("max_length"       , dict(type=int,
                                       default=100,
                                       metavar="N",
                                       help="Maximum read length permitted to be included"+
                                            " (BAM & bowtie files only. Default: %(default)s)")),
            ]

        big_genome = [
            ("big_genome"       , dict(action="store_true",
                                       default=False,
                                       help="Use slower but memory-efficient implementation "+
                                            "for big genomes or for memory-limited computers. For wiggle & bowtie files only.")),
            ]

        maxmem = [
            ("maxmem"           , dict(type=float,default=0,
                                       help="Maximum desired memory footprint in MB to devote to BigBed/BigWig files. May be exceeded by large queries. (Default: 0, No maximum)")), 
            ]
        
        # filetype-specific options
        self.filetype_options = {
             "BAM"    : length_ops,
             "bowtie" : length_ops + big_genome,
             "wiggle" : big_genome,
             "bigwig" : maxmem,
                                 
             }

        if self.allow_mapping == False:
            self.map_rules = []
            self.map_ops = []
        else:
            map_rules = [
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
                ]
            map_ops = [
                ("offset"           , dict(default=0,
                                            metavar="OFFSET",
                                            help="For `--fiveprime` or `--threeprime`, provide an integer "+
                                              "representing the offset into the read, "+
                                              "starting from either the 5\' or 3\' end, at which data "+
                                              "should be plotted. For `--fiveprime_variable`, "+
                                              "provide the filename of a two-column tab-delimited text "+
                                              "file, in which first column represents read length or the "+
                                              "special keyword `'default'`, and the second column represents "+
                                              "the offset from the five prime end of that read length at which the read should be mapped. (Default: %(default)s)")),
                ("nibble"           , dict(type=int,
                                            default=0,
                                            metavar="N",
                                            help="For use with `--center` only. nt to remove from each "+
                                                 "end of read before mapping (Default: %(default)s)")),
                ]
            
            for epoint in pkg_resources.iter_entry_points(group="plastid.mapping_rules"):
                reg_name = epoint.name
                pdict = epoint.load()
                if "name" in pdict:
                    reg_name = pdict.pop("name")
                    
                pdict["const"]  = reg_name
                pdict["action"] = "store_const"
                pdict["dest"]   = "%smapping" % self.prefix
                
                bamfunc = pdict.get("bamfunc",None)
                bowtiefunc = pdict.get("bowtiefunc",None)
                    
                if bamfunc is not None:
                    self.bamfuncs[reg_name] = bamfunc
                    pdict.pop("bamfunc")
                    
                if bowtiefunc is not None:
                    self.bowtiefuncs[reg_name] = bowtiefunc
                    pdict.pop("bowtiefunc")
                        
                map_rules.append((reg_name,pdict))
            
            for epoint in pkg_resources.iter_entry_points(group="plastid.mapping_options"):
                reg_name = epoint.name
                pdict = epoint.load()
                if "name" in pdict:
                    reg_name = pdict.pop("name")
                    
                map_ops.append((reg_name,pdict))
            
            self.map_rules = map_rules
            self.map_ops = map_ops
    
    def get_parser(self,
                   title=_DEFAULT_ALIGNMENT_FILE_PARSER_TITLE,
                   description=_DEFAULT_ALIGNMENT_FILE_PARSER_DESCRIPTION,
                   **kwargs):
        """Return an :py:class:`~argparse.ArgumentParser` that opens
        alignment (`BAM`_,  or `bowtie`_) or count (`Wiggle`_, `bedGraph`_) files.
         
        In the case of `bowtie`_ or `BAM`_ import, also parse arguments for mapping
        rules (e.g. fiveprime end mapping, threeprime end mapping, et c) and optional 
        read length filters
        
        
        Parameters
        ----------
        groupname : str, optional
            Name of argument group. If not `None`, an argument group with
            the specified name will be created and added to the parser.
            If not, arguments will be in the main group.         
        
        title : str, optional
            title for option group (used in command-line help screen)
                
        description : str, optional
            description of parser (used in command-line help screen)

        kwargs : keyword arguments
            Additional arguments to pass to :meth:`Parser.get_parser`

            
        Returns
        -------
        :class:`argparse.ArgumentParser`
        """        
        parser = Parser.get_parser(self,title=title,description=description,**kwargs)
        extra_args = []
        for k in self.input_choices:
            arglist = self.filetype_options.get(k,[])
            for arg in arglist:
                if arg not in extra_args:
                    extra_args.append(arg)
                
        if len(extra_args) > 0:
            # use mutator function- add new parser to `parser`
            Parser.get_parser(self,
                              parser=parser,
                              arglist=extra_args,title=title,description=description)
        
        if self.allow_mapping == True:
            Parser.get_parser(self,
                              parser=parser,
                              groupname="mapping_options",
                              arglist=self.map_rules,
                              title=_MAPPING_RULE_TITLE,
                              description=_MAPPING_RULE_DESCRIPTION,
                              )
            Parser.get_parser(self,
                              parser=parser,
                              groupname="sub_options",
                              arglist=self.map_ops,
                              title=_MAPPING_OPTION_TITLE,
                              description=_MAPPING_OPTION_DESCRIPTION,
                              )
        
        return parser
    
    def get_genome_array_from_args(self,args,printer=None):
        """Return a |GenomeArray|, |SparseGenomeArray| or |BAMGenomeArray|
        from arguments parsed by :py:func:`get_alignment_file_parser`
        
        Parameters
        ----------
        args : :py:class:`argparse.Namespace`
            Arguments from the parser
    
        printer : file-like, optional
            A stream to which stderr-like info can be written (default: |NullWriter|) 
        
        
        Returns
        -------
        |GenomeArray|, |SparseGenomeArray|, or |BAMGenomeArray|
        """
        from plastid.genomics.genome_array import GenomeArray, SparseGenomeArray,\
                                                   BAMGenomeArray, BigWigGenomeArray,\
                                                   SizeFilterFactory, CenterMapFactory,\
                                                   FivePrimeMapFactory, ThreePrimeMapFactory,\
                                                   VariableFivePrimeMapFactory,\
                                                   five_prime_map,  \
                                                   three_prime_map, \
                                                   center_map,      \
                                                   variable_five_prime_map        
        
        args = PrefixNamespaceWrapper(args,self.prefix)
        disabled = self.disabled
        map_rule = args.mapping
    
        if printer is None:
            printer = NullWriter()
        
        # require at least one countfile
        if len(args.count_files) == 0:
            printer.write("Please include at least one input file.")
            sys.exit(1)
        
        # require mapping rules unless wiggle
        if map_rule is None and args.countfile_format in ("BAM","bowtie"):
            printer.write("Please specify a read mapping rule.")
            sys.exit(1)
        
        if "countfile_format" not in disabled:
            
            if args.countfile_format in ("BAM","CRAM"):
                count_files = [pysam.Samfile(X,"rb") for X in args.count_files]
                try:
                    ga = BAMGenomeArray(count_files)
                except ValueError:
                    printer.write("Input BAM file(s) not indexed. Please index via:")
                    printer.write("")
                    for fn in args.count_files:
                        printer.write("    samtools index [-b|-c] %s" % fn)
                    printer.write("")
                    printer.write("Exiting.")
                    sys.exit(1)
                    
                size_filter = SizeFilterFactory(min=args.min_length,max=args.max_length)
                ga.add_filter("size:%s-%s" % (args.min_length,args.max_length) ,size_filter)
                if map_rule == "fiveprime":
                    map_function = FivePrimeMapFactory(int(args.offset))
                elif map_rule == "threeprime":
                    map_function = ThreePrimeMapFactory(int(args.offset))
                elif map_rule == "center":
                    map_function = CenterMapFactory(args.nibble)
                elif map_rule == "fiveprime_variable":
                    if str(args.offset) == "0":
                        printer.write("Please specify a filename to use for fiveprime variable offsets in --offset.")
                        sys.exit(1)
                    offset_dict = _parse_variable_offset_file(CommentReader(open(args.offset)))
                    map_function = VariableFivePrimeMapFactory(offset_dict)
                elif map_rule in self.bamfuncs:
                    map_function = functools.partial(self.bamfuncs[map_rule],args=args)
                else:
                    printer.write("Mapping rule '%s' not implemented for BAM input. Exiting." % map_rule)
                    sys.exit(1)
                ga.set_mapping(map_function)
                
            elif args.countfile_format == "bigwig":
                ga = BigWigGenomeArray(maxmem=args.maxmem)
                for align_file in args.count_files:
                    ga.add_from_bigwig("%s_fw.bw" % align_file,"+")
                    ga.add_from_bigwig("%s_rc.bw" % align_file,"-")
            
            # wiggle/bedGraph and bowtie
            else:
                if "big_genome" not in disabled and args.big_genome == True:
                    ga = SparseGenomeArray()
                else:
                    ga = GenomeArray()

                # wiggle/bedGraph
                if args.countfile_format == "wiggle":
                    for align_file in args.count_files:
                        printer.write("Opening wiggle files %s ..." % align_file)
                        with open("%s_fw.wig" % align_file) as fh:
                            ga.add_from_wiggle(fh,"+")
                        with open("%s_rc.wig" % align_file) as fh:
                            ga.add_from_wiggle(fh,"-")

                # bowtie
                elif args.countfile_format == "bowtie":
                    trans_args = { "nibble" : int(args.nibble) }
                    if map_rule == "fiveprime_variable":
                        transformation = variable_five_prime_map
                        if str(args.offset) == "0":
                            printer.write("Please specify a filename to use for fiveprime variable offsets in --offset.")
                            sys.exit(1)
                        else:
                            with open(args.offset) as myfile:
                                trans_args["offset"] = _parse_variable_offset_file(CommentReader(myfile))
                    else:
                        trans_args["offset"] = int(args.offset)
                        if map_rule == "fiveprime":
                            transformation = five_prime_map
                        elif map_rule == "threeprime":
                            transformation = three_prime_map
                        elif map_rule == "entire":
                            transformation = center_map
                        elif map_rule == "center":
                            transformation = center_map
                        elif  map_rule in self.bowtiefuncs:
                            transformation = self.bowtiefuncs[map_rule]
                            trans_args["args"] = args
                        else:
                            printer.write("Mapping rule '%s' not implemented for bowtie input. Exiting." % map_rule)
                            sys.exit(1)
                
                    for infile in args.count_files:
                        with opener(infile) as my_file:
                            ga.add_from_bowtie(my_file,transformation,min_length=args.min_length,max_length=args.max_length,**trans_args)
                
        printer.write("Counted %s total reads." % ga.sum())
        
        if "sum" not in disabled and args.sum is not None:
            ga.set_sum(args.sum)
            
        if "normalize" not in disabled and args.normalize == True:
            printer.write("Normalizing to reads per million.")
            ga.set_normalize(True)
        
        return ga


#===============================================================================
# INDEX: Annotation file parser
#===============================================================================

class AnnotationParser(Parser):
    """Parser for annotation files in various formats
    
    Parameters
    ----------
    groupname : str, optional
        Name of argument group. If not `None`, an argument group with
        the specified name will be created and added to the parser.
        If not, arguments will be in the main group.         
    
    prefix : str, optional
        string prefix to add to default argument options (Default: "")

    disabled : list, optional
        list of parameter names that should be disabled from parser,
        without preceding dashes

    input_choices : list, optional
        list of permitted alignment file type choices for input
    
    allow_mapping : bool, optional
        Enable/disable user configuration of mapping rules (default: True)    
    """
    
    def __init__(self,
                 prefix="",
                 disabled=None,
                 groupname="annotation_options",
                 input_choices=("BED","BigBed","GTF2","GFF3")
                ):
        """Create a parser for genomic features in an annotation file
        
        Parameters
        ----------
        groupname : str, optional
            Name of argument group. If not `None`, an argument group with
            the specified name will be created and added to the parser.
            If not, arguments will be in the main group.         
        
        prefix : str, optional
            string prefix to add to default argument options (Default: "")

        disabled : list, optional
            list of parameter names that should be disabled from parser,
            without preceding dashes

        input_choices : list, optional
            list of permitted alignment file type choices for input
        
        allow_mapping : bool, optional
            Enable/disable user configuration of mapping rules (default: True)
        """
        Parser.__init__(self,groupname=groupname,prefix=prefix,disabled=disabled)
        self.input_choices = input_choices
        self.arguments = [
                ("annotation_files"     , dict(metavar="infile.[%s]" % " | ".join(input_choices),# | psl]",
                                              type=str,nargs="+",default=[],
                                              help="Zero or more annotation files (max 1 file if BigBed)")),                               
                ("annotation_format"   , dict(choices=input_choices,
                                              default="GTF2",
                                              help="Format of %sannotation_files (Default: GTF2). Note: GFF3 assembly assumes SO v.2.5.2 feature ontologies, which may or may not match your specific file." % prefix)),    
                ("add_three"           , dict(default=False,
                                              action="store_true",
                                              help="If supplied, coding regions will be extended by 3 nucleotides at their 3\' ends (except for GTF2 files that explicitly include `stop_codon` features). Use if your annotation file excludes stop codons from CDS.")),
                ("tabix"               , dict(default=False,
                                              action="store_true",
                                              help="%sannotation_files are tabix-compressed and indexed (Default: False). Ignored for BigBed files." % prefix)),
                ("sorted"              , dict(default=False,
                                              action="store_true",
                                              help="%sannotation_files are sorted by chromosomal position (Default: False)" % prefix)),
            ]
        
        # options for specific filetypes
        self.filetype_options = {
            "BED" : [("bed_extra_columns", dict(default=0,
                                                nargs="+",
                                                help="Number of extra columns in BED file (e.g. in custom ENCODE formats) "+
                                                     "or list of names for those columns. (Default: %(default)s)."))
                    ],
            "BigBed" : [
                ("maxmem"           , dict(type=float,default=0,
                                           help="Maximum desired memory footprint in MB to devote to BigBed/BigWig files. May be exceeded by large queries. (Default: 0, No maximum)")), 
                        
                        ],
            "GFF3" : [("gff_transcript_types", dict(type=str,
                                                    default=_DEFAULT_GFF3_TRANSCRIPT_TYPES,
                                                    nargs="+",
                                                    help="GFF3 feature types to include as transcripts, even "+\
                                                         "if no exons are present (for GFF3 only; default: use SO v2.5.3 specification)")),
                      ("gff_exon_types", dict(type=str,
                                              default=_DEFAULT_GFF3_EXON_TYPES,
                                              nargs="+",
                                              help="GFF3 feature types to include as exons (for GFF3 only; default: use SO v2.5.3 specification)")),
                      ("gff_cds_types", dict(type=str,
                                             default=_DEFAULT_GFF3_CDS_TYPES,
                                             nargs="+",
                                             help="GFF3 feature types to include as CDS (for GFF3 only; default: use SO v2.5.3 specification)")),
                     ]
             }

    def get_parser(self,
                   title=_DEFAULT_ANNOTATION_PARSER_TITLE,
                   description=_DEFAULT_ANNOTATION_PARSER_DESCRIPTION,
                   **kwargs):
        """Return an :class:`~argparse.ArgumentParser` that opens annotation files.
         
        Parameters
        ----------
        title : str, optional
            title for option group (used in command-line help screen)
                
        description : str, optional
            description of parser (used in command-line help screen)

        kwargs : keyword arguments
            Additional arguments to pass to :meth:`Parser.get_parser`
            
        Returns
        -------
        :class:`argparse.ArgumentParser`
        """        
        parser = Parser.get_parser(self,title=title,description=description,**kwargs)
        
        for k in self.input_choices:
            arglist = self.filetype_options.get(k)
            if arglist is not None:
                # use mutator function- add new parser to `parser`
                Parser.get_parser(self,
                                  parser=parser,
                                  groupname="%s_%s_options" % (self.groupname,k),
                                  title="%s-specific options" % k,arglist=arglist)
        
        return parser           

    def get_transcripts_from_args(self,args,printer=None,return_type=None,require_sort=False):
        """Return a generator of |Transcript| objects from arguments parsed by :func:`get_annotation_file_parser`
        
        Parameters
        ----------
        args : :py:class:`argparse.Namespace`
            Namespace object from :py:func:`get_annotation_file_parser`
        
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

        return self.get_segmentchains_from_args(args, printer=printer, return_type=Transcript, require_sort=require_sort)

    def get_segmentchains_from_args(self,args,printer=None,return_type=None,require_sort=False):
        """Return a generator of |SegmentChain| objects from arguments parsed by :py:func:`get_annotation_file_parser`
        
        Parameters
        ----------
        args : :py:class:`argparse.Namespace`
            Namespace object from :py:func:`get_annotation_file_parser`
        
        printer : file-like, optional
            A stream to which stderr-like info can be written (Default: |NullWriter|) 
        
        return_type : |SegmentChain| or subclass, optional
            Type of object to return (Default: |Transcript|)
    
        require_sort : bool, optional
            If True, quit if the annotation file(s) are not sorted or indexed
        
        Returns
        -------
        iterator
            |SegmentChain| objects, either in order of appearance (if input was a
            `BED`_, `BigBed`_, or `PSL`_ file), or sorted lexically by chromosome,
            start coordinate, end coordinate, and then strand (if input was `GTF2`_
            or `GFF3`_).
        
        
        See Also
        --------
        get_annotation_file_parser
            Function that creates :py:class:`argparse.ArgumentParser` whose output
            :py:class:`~argparse.Namespace` is processed by this function    
        """
        if printer is None:
            printer = NullWriter()
            
        if return_type is None:
            return_type = SegmentChain
    
        args = PrefixNamespaceWrapper(args,self.prefix)
        disabled = self.disabled
    
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
                    printer.write(GFF_SORT_MESSAGE.replace("FORMAT",args.annotation_format))
                    sys.exit(1)
    
        printer.write("Parsing features in %s ..." % ", ".join(args.annotation_files))
        
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
            transcripts = iter(BigBedReader(args.annotation_files[0],
                                       return_type=return_type,
                                       cache_depth=1,
                                       add_three_for_stop=add_three,
                                       printer=printer,
                                       maxmem=args.maxmem))
            
        elif tabix == True:
            #streams = [pysam.tabix_iterator(opener(X), lambda x,y: x) for X in args.annotation_files] # used to work in earlier pysam
            # string parsing by supplying None instead of `asTuple()` no longer works
            # nor do anonymous lambda functions
            streams = [pysam.tabix_iterator(opener(X), pysam.asTuple()) for X in args.annotation_files]
        else:
            streams = (opener(X) for X in args.annotation_files)
    
        if args.annotation_format in ("GFF3","GTF2"):
            from plastid.readers.gff import GFF3_TranscriptAssembler, GTF2_TranscriptAssembler
            if 'sorted' not in disabled and args.sorted == False and \
               'tabix' not in disabled and args.tabix == False and \
               any((os.stat(X).st_size >= _GFF_SORT_SIZE for X in args.annotation_files)):
                msg = """Transcript assembly on large FORMAT files can require a lot of memory.
    Consider using a sorted file with the '--sorted' flag and/or tabix-compression.
    """
                msg += GFF_SORT_MESSAGE
                msg = msg.replace("FORMAT",args.annotation_format)
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
        
    def get_genome_hash_from_args(self,args,printer=None):
        """Return a |GenomeHash| of regions from command-line arguments
    
        Parameters
        ----------
        args : :py:class:`argparse.Namespace`
            Namespace object from :py:func:`get_mask_file_parser`
        
        printer : file-like
            A stream to which stderr-like info can be written (Default: |NullWriter|) 
    
    
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
        from plastid.readers.bed import BED_Reader
        from plastid.readers.gff import GTF2_Reader, GFF3_Reader
        from plastid.readers.psl import PSL_Reader
        if printer is None:
            printer = NullWriter()
            
        prefix = self.prefix
        args = PrefixNamespaceWrapper(args,prefix)
        
        if len(args.annotation_files) > 0:
            printer.write("Opening mask annotation file(s) %s ..." % ", ".join(args.annotation_files))
            if args.annotation_format in ("BED","GTF2","GFF3") and args.tabix == False:
                msg = """Unindexed mask files can require lots of memory in large (e.g. mammalian) genomes.
    Consider converting to BigBed or using tabix to index your mask file."""
                warnings.warn(msg,ArgumentWarning)
    
            if len(args.annotation_files) > 0:
                if args.annotation_format == "BigBed":
                    if len(args.annotation_files) > 1:
                        printer.write("Bad arguments: we can only process one BigBed file.")
                        sys.exit(2)
                    return BigBedGenomeHash(args.annotation_files[0])
                elif "tabix" not in self.disabled and args.tabix == True:
                    return TabixGenomeHash(args.annotation_files,args.annotation_format,printer=printer)
                else:
                    streams = (opener(X) for X in args.annotation_files)
                    if args.annotation_format == "BED":
                        reader = BED_Reader
                    elif args.annotation_format == "GTF2":
                        reader = GTF2_Reader
                    elif args.annotation_format == "GFF3":
                        reader = GFF3_Reader
                    elif args.annotation_format == "PSL":
                        reader = PSL_Reader
                    else:
                        assert False
                        
                    hash_ivcs = list(reader(*streams))

                    return GenomeHash(hash_ivcs)
        else:
            return GenomeHash()



class MaskParser(AnnotationParser):
    """Create a parser for masking genomic features given in an annotation file
    
    Parameters
    ----------
    groupname : str, optional
        Name of argument group. If not `None`, an argument group with
        the specified name will be created and added to the parser.
        If not, arguments will be in the main group.         
    
    prefix : str, optional
        string prefix to add to default argument options (Default: "")

    disabled : list, optional
        list of parameter names that should be disabled from parser,
        without preceding dashes

    input_choices : list, optional
        list of permitted alignment file type choices for input
    
    allow_mapping : bool, optional
        Enable/disable user configuration of mapping rules (default: True)
    """

    def __init__(self,
                 prefix="mask_",
                 disabled=None,
                 groupname="mask_options",
                 input_choices=("BED","BigBed","GTF2","GFF3","PSL")
                ):
        """Create a parser for genomic features in an annotation file
        
        Parameters
        ----------
        groupname : str, optional
            Name of argument group. If not `None`, an argument group with
            the specified name will be created and added to the parser.
            If not, arguments will be in the main group.         
        
        prefix : str, optional
            string prefix to add to default argument options (Default: "")

        disabled : list, optional
            list of parameter names that should be disabled from parser,
            without preceding dashes

        input_choices : list, optional
            list of permitted alignment file type choices for input
        
        allow_mapping : bool, optional
            Enable/disable user configuration of mapping rules (default: True)
        """
        AnnotationParser.__init__(self,
                                  prefix=prefix,
                                  disabled=disabled,
                                  groupname=groupname,
                                  input_choices=("BED","BigBed","GTF2","GFF3","PSL"))

    def get_parser(self,
                   title=_MASK_PARSER_TITLE,
                   description=_MASK_PARSER_DESCRIPTION,
                   **kwargs):
        """Return an :py:class:`~argparse.ArgumentParser` that opens annotation files as masks
        alignment (`BAM`_ or `bowtie`_) or count (`Wiggle`_, `bedGraph`_) files.
         
        Parameters
        ----------
        title : str, optional
            title for option group (used in command-line help screen)
                
        description : str, optional
            description of parser (used in command-line help screen)

        arglist : list, optional
            If not `None`, arguments in this list will be added to `parser`.
            Otherwise, arguments will be taken from `self.arguments`.
            
            The list should be a list of tuples of ('argument_name',dict_of_options),
            where `argument_name` is a string, and `dict_of_options` a dictionary
            of keyword arguments to pass to :meth:`argparse.ArgumentParser.add_argument`.
            
        Returns
        -------
        :class:`argparse.ArgumentParser`
        """        
        return AnnotationParser.get_parser(self,title=title,description=description,**kwargs)
        
#===============================================================================
# INDEX: Sequence parser
#===============================================================================
        
class SequenceParser(AnnotationParser):
    """Parser for sequence files
            
    Parameters
    ----------
    groupname : str, optional
        Name of argument group. If not `None`, an argument group with
        the specified name will be created and added to the parser.
        If not, arguments will be in the main group.         
    
    prefix : str, optional
        string prefix to add to default argument options (Default: "")

    disabled : list, optional
        list of parameter names that should be disabled from parser,
        without preceding dashes

    input_choices : list, optional
        list of permitted alignment file type choices for input
    """
    
    def __init__(self,
                 groupname="sequence_options",
                 prefix="",
                 disabled=None,
                 input_choices=("fasta","fastq","twobit","genbank","embl"),
                 ):
        """Create a parser for genomic sequence
        
        Parameters
        ----------
        groupname : str, optional
            Name of argument group. If not `None`, an argument group with
            the specified name will be created and added to the parser.
            If not, arguments will be in the main group.         
        
        prefix : str, optional
            string prefix to add to default argument options (Default: "")

        disabled : list, optional
            list of parameter names that should be disabled from parser,
            without preceding dashes

        input_choices : list, optional
            list of permitted alignment file type choices for input
        """
        Parser.__init__(self,groupname=groupname,prefix=prefix,disabled=disabled)
        self.input_choices = input_choices
        self.arguments = [
                ("sequence_file"     , dict(metavar="infile.[%s]" % " | ".join(input_choices),
                                            type=str,
                                            help="A file of DNA sequence")),                               
                ("sequence_format"   , dict(choices=input_choices,
                                            default="fasta",
                                            help="Format of %ssequence_file (Default: fasta)." % prefix)),    
            ]
        
    def get_parser(self,
                   title=_DEFAULT_SEQUENCE_PARSER_TITLE,
                   description=_DEFAULT_SEQUENCE_PARSER_DESCRIPTION,
                   **kwargs):
        """Return an :py:class:`~argparse.ArgumentParser` that opens sequence files
         
        Parameters
        ----------
        
        title : str, optional
            title for option group (used in command-line help screen)
            
        description : str, optional
            description of parser (used in command-line help screen)

        kwargs : keyword arguments
            Additional arguments to pass to :meth:`Parser.get_parser`

            
        Returns
        -------
        :class:`argparse.ArgumentParser`
        
        
        See also
        --------
        get_seqdict_from_args
            function that parses the :py:class:`~argparse.Namespace` returned
            by this :py:class:`~argparse.ArgumentParser`
        """        
        
        return Parser.get_parser(self,title=title,description=description,**kwargs)

    def get_seqdict_from_args(self,args,index=True,printer=None):
        """Retrieve a dictionary-like object of sequences
    
        Parameters
        ----------
        args : :py:class:`argparse.Namespace`
            Namespace object from :py:func:`get_sequence_file_parser`
        
        index : bool, optional
            If sequence format is anything other than twobit, open with
            lazily-evaluating :func:`Bio.SeqIO.index` instead of
            :func:`Bio.SeqIO.to_dict` (Default: `True`)
            
        printer : file-like
            A stream to which stderr-like info can be written (Default: |NullWriter|) 
    
        Returns
        -------
        dict-like
            Dictionary-like object mapping chromosome names to
            :class:`Bio.SeqRecord.SeqRecord`-like objects
        """
        if printer is None:
            printer = NullWriter()
            
        args = PrefixNamespaceWrapper(args,self.prefix)
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
# INDEX: Plotting parser
#===============================================================================

class PlottingParser(Parser):
    """Parser for plotting options
        
    Parameters
    ----------
    groupname : str, optional
        Name of argument group. If not `None`, an argument group with
        the specified name will be created and added to the parser.
        If not, arguments will be in the main group.         
    
    prefix : str, optional
        string prefix to add to default argument options (Default: "")

    disabled : list, optional
        list of parameter names that should be disabled from parser,
        without preceding dashes    
    """

    def __init__(self,
                 groupname="plotting_options",
                 prefix="",
                 disabled=None):
        """Create a parser for plotting arguments
        
        Parameters
        ----------
        groupname : str, optional
            Name of argument group. If not `None`, an argument group with
            the specified name will be created and added to the parser.
            If not, arguments will be in the main group.         
        
        prefix : str, optional
            string prefix to add to default argument options (Default: "")

        disabled : list, optional
            list of parameter names that should be disabled from parser,
            without preceding dashes
        """        
        Parser.__init__(self,groupname=groupname,prefix=prefix,disabled=disabled)
        from matplotlib.backend_bases import FigureCanvasBase as fcb
        if len(prefix) > 0:
            prefix += "_"
    
        try:
            filetypes = sorted(fcb.get_supported_filetypes().keys())
            default_ftype = fcb.get_default_filetype()
        except: # matplotlib < 1.4.0
            filetypes = ["eps","jpeg","pdf","png","svg"]
            default_ftype = "pdf"
        
        self.arguments = [
                ("figformat",dict(default=default_ftype,type=str,choices=filetypes,
                                  help="File format for figure(s); Default: %(default)s)")),
                ("figsize",  dict(nargs=2,default=None,type=float,metavar="N",
                                  help="Figure width and height, in inches. (Default: use matplotlibrc params)"
                                  )),
                ("title",    dict(type=str,default=None,help="Base title for plot(s).")),
                ("cmap",     dict(type=str,default=None,
                                  help="Matplotlib color map from which palette will be made (e.g. 'Blues','autumn','Set1'; default: use colors from ``--stylesheet`` if given, or color cycle in matplotlibrc)"
                                  )),
                ("dpi",      dict(type=int,default=150,
                                  help="Figure resolution (Default: %(default)s)")),
            ]
        
        try:
            import matplotlib.style
            stylesheets = matplotlib.style.available
    
            if "stylesheet" not in self.disabled:
                self.arguments.append(("stylesheet",
                                       dict(default=None,choices=stylesheets,
                                            help="Use this matplotlib stylesheet instead of matplotlibrc params")
                                     ))
        except ImportError: # matplotlib < 1.4.0
            pass
            
    def get_parser(self,title=_DEFAULT_PLOTTING_TITLE,description=None):
        """Return an :py:class:`~argparse.ArgumentParser` to control plotting     
    
        Parameters
        ----------
            
        title : str, optional
            title for option group (used in command-line help screen)
            
        description : str, optional
            description of parser (used in command-line help screen)
        
       
        Returns
        -------
        :class:`argparse.ArgumentParser`
        """
        return Parser.get_parser(self,title=title,description=description)


    def set_style_from_args(self,args):
        """Parse style information, if present on system and defined in `args`

        Parameters
        ----------
        args : :class:`argparse.Namespace`
            Namespace object from :func:`get_plotting_parser`
        """
        try:
            import matplotlib.style
            if getattr(args,"stylesheet",None) is not None:
                matplotlib.style.use(args.stylesheet)
        except ImportError:
            pass

    def get_figure_from_args(self,args,**kwargs):
        """Return a :class:`matplotlib.figure.Figure` following arguments from :func:`get_plotting_parser`
    
        A new figure is created with parameters specified in `args`. If these are 
        not found, values found in `**kwargs` will instead be used. If these are 
        not found, we fall back to matplotlibrc values.
    
        Parameters
        ----------
        args : :class:`argparse.Namespace`
            Namespace object from :func:`get_plotting_parser`
    
        kwargs : keyword arguments
            Fallback arguments for items not defined in `args`, plus any other
            keyword arguments.
    
        Returns
        -------
        :class:`matplotlib.figure.Figure`
            Matplotlib figure
        """
        import matplotlib.pyplot as plt
        args = PrefixNamespaceWrapper(args,self.prefix)
    
        fargs = {}
        # keep this loop in place in case we add additional command line attributes as fig properties later
        for attr in ("figsize",): #,"dpi"): # dpi if applied in plt.figure() doesn't get used in saving; 
            if attr in kwargs:
                v = kwargs[attr]
            else:
                v = getattr(args,attr,None)

            if v is not None:
                fargs[attr] = v
    
        # copy values from fargs
        kwargs.update(fargs)
        return plt.figure(**kwargs)

    def get_colors_from_args(self,args,num_colors):
        """Return a list of colors from arguments parsed by a parser from :func:`get_plotting_parser`
    
        If a matplotlib colormap is specified in `args.figcolors`, colors will be
        generated from that map.
    
        Otherwise, if a stylesheet is specified, colors will be fetched from 
        the stylesheet's color cycle.
        
        Otherwise, colors will be chosen from the default color cycle specified
        ``matplotlibrc``.
    
    
        Parameters
        ----------
        args : :class:`argparse.Namespace`
            Namespace object from :func:`get_plotting_parser`
    
        num_colors : int
            Number of colors to fetch
        
    
        Returns
        -------
        list
            List of matplotlib colors
        """
        import matplotlib.cm
        args       = PrefixNamespaceWrapper(args,self.prefix)
        figcolors  = getattr(args,"cmap",None)
    
        if figcolors is not None:
            import numpy
            cmap = matplotlib.cm.get_cmap(figcolors) 
            if num_colors > 1:
                colors = cmap(numpy.linspace(0,1.0,num_colors))
            else:
                colors = [cmap(0.5)]
        else:   
            from itertools import cycle
            try:
                
                color_cycle = cycle(matplotlib.rcParams["axes.prop_cycle"].by_key()["color"])
            except KeyError:
                color_cycle = cycle(matplotlib.rcParams["axes.color_cycle"])

            colors = [next(color_cycle) for _ in range(num_colors)]
    
        return colors


#===============================================================================
# INDEX: Parser for generic command-line options (e.g. warning control)
#===============================================================================

class BaseParser(Parser):
    """Parser basic options
    
    Parameters
    ----------
    groupname : str, optional
        Name of argument group. If not `None`, an argument group with
        the specified name will be created and added to the parser.
        If not, arguments will be in the main group.         
    
    prefix : str, optional
        string prefix to add to default argument options (Default: "")

    disabled : list, optional
        list of parameter names that should be disabled from parser,
        without preceding dashes    
    """
    
    def __init__(self,
                 groupname="base_options",
                 prefix="",
                 disabled=None,
                 ):
        """Create a parser for basic options for command-line scripts, such as warnings and logging
        
        Parameters
        ----------
        groupname : str, optional
            Name of argument group. If not `None`, an argument group with
            the specified name will be created and added to the parser.
            If not, arguments will be in the main group.         
        
        prefix : str, optional
            string prefix to add to default argument options (Default: "")

        disabled : list, optional
            list of parameter names that should be disabled from parser,
            without preceding dashes
        """
        Parser.__init__(self,groupname=groupname,prefix=prefix,disabled=disabled)
        self.arguments = []
#         self.level_desc = ["--silent","--quiet","--verbose","--raise"]
        
    def get_parser(self,title=None,description=None):
        """Return an :py:class:`~argparse.ArgumentParser`     
    
        Parameters
        ----------
            
        title : str, optional
            title for option group (used in command-line help screen)
            
        description : str, optional
            description of parser (used in command-line help screen)
        
       
        Returns
        -------
        :class:`argparse.ArgumentParser`
        """
        p = Parser.get_parser(self)
        g = p.add_argument_group(title="warning/error options")
        
        
        g.add_argument("-q","--quiet",dest="warnlevel",action="store_const",const=-1,
                       help="Suppress all warning messages. Cannot use with '-v'.")
        g.add_argument("-v","--verbose",dest="warnlevel",action="count",
                       help="Increase verbosity. With '-v', show every warning. With '-vv', turn warnings into exceptions. Cannot use with '-q'. (Default: show each type of warning once)")
        
#         g.add_argument("--silent",dest="warnlevel",action="store_const",const=0,
#                        help="Suppress all warning messages")
#         g.add_argument("--quiet",dest="warnlevel",action="store_const",const=1,
#                        help="Show each type of warning once (Default)")
#         g.add_argument("--verbose",dest="warnlevel",action="store_const",const=2,
#                        help="Show every warning instance")
#         g.add_argument("--raise",dest="warnlevel",action="store_const",const=3,
#                        help="Raise exceptions instead of warnings")
        p.set_defaults(warnlevel=0)

        return p

    def get_base_ops_from_args(self,args):
        
        global warnings
        
        args = PrefixNamespaceWrapper(args,self.prefix)
        warnlevel = args.warnlevel
        actions = ["ignore",
                   "onceperfamily",
                   "always",
                   "error"]
#         desc = self.level_desc
        
#         if len(set(desc) & set(sys.argv[1:])) > 1:
#             warnings.warn("Only one of [%s] is permitted. Using %s" % (", ".join(desc),desc[warnlevel]),
#                            ArgumentWarning)
            
        if warnlevel >= len(actions) - 1:
            warnlevel = len(actions) - 2
        try:
            action = actions[warnlevel+1]
        except IndexError:
            warnings.warn("Invalid warning level. Expected 0-3, found %s. Defaulting to level 1 (`--once`)." % warnlevel,UserWarning)
            action = actions[1]
        
        for type_, msg in PLASTID_WARNINGS:
            filterwarnings(action,message=msg,category=type_)


PLASTID_WARNINGS = [
                    
    # mapping rules
    (DataWarning,"File contains read alignments shorter"),
    (DataWarning,"No offset for reads of length"),
    (DataWarning,"longer than read length"),
    
    # genome_array
    (DataWarning,"Temporarily turning off normalization"),
    
    # roi_tools
    (DataWarning,"is a zero-length SegmentChain. Returning 0-length count vector"),
    
    # metagene
    (Warning,r"IndexError finding common positions at region.*"),
    (DataWarning,"has no gene_id. Inferring gene_id"),
    (DataWarning,"has no attribute"),
    (DataWarning,"Ignoring labels"),
    
    # phase_by_size
    (DataWarning,"is not divisible by 3. Ignoring last partial codon."),
    
    # util.io.filters
    (Warning,"Could not alert listener"),

    # util.services.decorators
    (DeprecationWarning,"is deprecated and will be removed from module"),
    
    # gff
    (DataWarning,"because it contains exons on multiple chromosomes or strands"),
    (DataWarning,"because start or stop codons are outside exon boundaries"),
    (DataWarning,"with no `Parent` or `ID`. Ignoring."),
    (DataWarning,"because it contains exons on multiple strands"),
    (DataWarning,"because start or stop codons are outside exon boundaries."),

    # bed
    (FileFormatWarning,"Extra columns specified by."),
    (FileFormatWarning,"Are you sure this is a"),
    (FileFormatWarning,"Are you sure this BED file has extra columns"),
    (FileFormatWarning,"Maybe this BED has extra columns"),
    
    # gff_tokens
    (FileFormatWarning,"Found duplicate attribute key"),

    # BigBed
    (FileFormatWarning,"Could not find or could not parse autoSql declaration in BigBedFile"),
    
    # autoSql
    (DataWarning,"Could not convert autoSql value"),
    
    # psl
    (FileFormatWarning,"Rejecting line")

]


#===============================================================================
# INDEX: Deprecated alignment functions, now aliased to classes above
#===============================================================================


@deprecated(version="0.5.0",instead="AlignmentParser")
def get_alignment_file_parser(input_choices=("BAM","bigwig","bowtie","wiggle"),
                              disabled=None,
                              prefix="",
                              title=_DEFAULT_ALIGNMENT_FILE_PARSER_TITLE,
                              description=_DEFAULT_ALIGNMENT_FILE_PARSER_DESCRIPTION,
                              map_desc=_MAPPING_RULE_DESCRIPTION,
                              return_subparsers=False):
    tmp = AlignmentParser(input_choices=input_choices,prefix=prefix,disabled=disabled)
    return tmp.get_parser(title=title, description=description)

@deprecated(version="0.5.0",instead="AlignmentParser.get_genome_array_from_args()")
def get_genome_array_from_args(args,prefix="",disabled=None,printer=None):
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
    tmp = AlignmentParser(prefix=prefix,disabled=disabled)
    return tmp.get_genome_array_from_args(args,printer=printer)


#===============================================================================
# INDEX: deprecated annotation file parser, and helper functions
#===============================================================================

@deprecated(version="0.5.0",instead="AnnotationParser")
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
    tmp = AnnotationParser(groupname="annotation_options",
                           prefix=prefix,
                           disabled=disabled,
                           input_choices=input_choices)
    parser = tmp.get_parser(title, description)
    return parser

@deprecated(version="0.5.0",instead="AnnotationParser.get_transcripts_from_args()")
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
    tmp = AnnotationParser(groupname="annotation_options",
                       prefix=prefix,
                       disabled=disabled)
    return tmp.get_transcripts_from_args(args,
                                         printer=printer,
                                         return_type=return_type,
                                         require_sort=require_sort)

@deprecated(version="0.5.0",instead="AnnotationParser.get_parser()")
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

@deprecated(version="0.5.0",instead="AnnotationParser.get_transcripts_from_args()")
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
        A stream to which stderr-like info can be written (Default: |NullWriter|) 
    
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
    disabled.append([prefix+"add_three"])
    return get_transcripts_from_args(args,
                                     prefix=prefix,
                                     disabled=disabled,
                                     printer=printer,
                                     return_type=SegmentChain,
                                     require_sort=require_sort)

@deprecated(version="0.5.0",instead="AnnotationParser")
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
    tmp = AnnotationParser(groupname="%s_options" % prefix,
                           prefix=prefix,
                           disabled=disabled,
                           input_choices=["BED","GTF2","GFF3","BigBed","PSL"])
    return tmp.get_parser(_MASK_PARSER_TITLE, _MASK_PARSER_DESCRIPTION) 

@deprecated(version="0.5.0",instead="AnnotationParser.get_genome_hash()")
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
        A stream to which stderr-like info can be written (Default: |NullWriter|) 


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
    tmp = AnnotationParser(groupname="mask_options",prefix=prefix)
    return tmp.get_genome_hash_from_args(args, printer=printer)

#===============================================================================
# INDEX: deprecated sequence file parser
#===============================================================================

@deprecated(version="0.5.0",instead="SequenceParser")
def get_sequence_file_parser(input_choices=("fasta","fastq","twobit","genbank","embl"),
                             disabled=(),
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
    tmp = SequenceParser(disabled=disabled,prefix=prefix,input_choices=input_choices)
    return tmp.get_parser(title=title,description=description)

@deprecated(version="0.5.0",instead="SequenceParser.get_seqdict_from_args()")
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
        A stream to which stderr-like info can be written (Default: |NullWriter|) 

    Returns
    -------
    dict-like
        Dictionary-like object mapping chromosome names to
        :class:`Bio.SeqRecord.SeqRecord`-like objects
    """
    tmp = SequenceParser(prefix=prefix)
    return tmp.get_seqdict_from_args(args, index=index, printer=printer)


#===============================================================================
# INDEX: deprecated plotting
#===============================================================================


@deprecated(version="0.5.0",instead="PlottingParser")
def get_plotting_parser(prefix="",disabled=[],title=_DEFAULT_PLOTTING_TITLE):
    """Return an :py:class:`~argparse.ArgumentParser` to control plotting     

    Parameters
    ----------
        
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
    get_colors_from_args
        parse colors and/or colormaps from this argument parser
    """
    tmp = PlottingParser(prefix=prefix,disabled=disabled)
    return tmp.get_parser(title=title)

@deprecated(version="0.5.0",instead="PlottingParser.get_figure_from_args()")
def get_figure_from_args(args,**kwargs):
    """Return a :class:`matplotlib.figure.Figure` following arguments from :func:`get_plotting_parser`

    A new figure is created with parameters specified in `args`. If these are 
    not found, values found in `**kwargs` will instead be used. If these are 
    not found, we fall back to matplotlibrc values.

    Parameters
    ----------
    args : :class:`argparse.Namespace`
        Namespace object from :func:`get_plotting_parser`

    kwargs : keyword arguments
        Fallback arguments for items not defined in `args`, plus any other
        keyword arguments.

    Returns
    -------
    :class:`matplotlib.figure.Figure`
        Matplotlib figure
    """
    tmp = PlottingParser()
    return tmp.get_figure_from_args(args,**kwargs)

@deprecated(version="0.5.0",instead="PlottingParser.get_colors_from_args()")
def get_colors_from_args(args,num_colors):
    """Return a list of colors from arguments parsed by a parser from :func:`get_plotting_parser`

    If a matplotlib colormap is specified in `args.figcolors`, colors will be
    generated from that map.

    Otherwise, if a stylesheet is specified, colors will be fetched from 
    the stylesheet's color cycle.
    
    Otherwise, colors will be chosen from the default color cycle specified
    ``matplotlibrc``.


    Parameters
    ----------
    args : :class:`argparse.Namespace`
        Namespace object from :func:`get_plotting_parser`

    num_colors : int
        Number of colors to fetch
    

    Returns
    -------
    list
        List of matplotlib colors
    """
    tmp = PlottingParser()
    return tmp.get_colors_from_args(args,num_colors)


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
