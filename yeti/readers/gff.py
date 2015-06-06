#!/usr/bin/env python
"""Tools for reading, writing, analyzing, and manipulating GFFs and features.

Convenience functions
---------------------
:py:func:`GTF2_to_Transcripts`
    Convenience method to assemble entries in GTF2 files into |Transcript|

:py:func:`GFF3_to_Transcripts`
    Convenience method to assemble entries in GFF3 files into |Transcript|


Important classes
-----------------
|GTF2_Reader|
    Read single |SegmentChain| objects line-by-line from GTF2 files
    
|GFF3_Reader|
    Read single |SegmentChain| objects line-by-line from GFF3 files

|GTF2_TranscriptAssembler|
    Assembles |Transcript| objects from one or more features in one or more GTF2 files

|GFF3_TranscriptAssembler|
    Assembles |Transcript| objects from one or more features in one or more GFF3 files


See Also
--------
`GFF3 specification <http://song.sourceforge.net/gff3.shtml>`_
    GFF3 specification by the Sequence Ontology consortium

`GTF2.2 specification <http://mblab.wustl.edu/GTF22.html>`_
    Hosted by the Brent lab
    
`UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_.
    GFF & GTF descriptions at UCSC
"""
__author__="joshua"
__date__ ="$Dec 1, 2010 11:00:55 AM$"
import datetime
import itertools
import gc
import copy
from abc import abstractmethod
from yeti.util.io.filters import AbstractWriter, AbstractReader, SkipBlankReader
from yeti.util.io.openers import NullWriter
from yeti.readers.common import add_three_for_stop_codon, \
                                                get_identical_attributes, \
                                                AssembledFeatureReader
from yeti.genomics.roitools import Transcript, SegmentChain, \
                                         GenomicSegment, \
                                         sort_segmentchains_lexically
from yeti.util.services.decorators import deprecated
from yeti.readers.gff_tokens import _make_generic_tokens, \
                                                    make_GFF3_tokens, \
                                                    make_GTF2_tokens, \
                                                    parse_GFF3_tokens, \
                                                    parse_GTF2_tokens


#===============================================================================
# INDEX: SO v2.5.3 feature types
#   see: http://www.sequenceontology.org/resources/intro.html
#===============================================================================

_DEFAULT_GFF3_GENE_TYPES = {
    "gene",
    "candidate_gene",
    "functional_candidate_gene",
    "positional_candidate_gene",
    "cryptic_gene",
    "cryptogene",
    "engineered_gene",
    "engineered_foreign_gene",
    "engineered_foreign_transposable_element_gene",
    "engineered_fusion_gene"
    "epigenetically_modified_gene",
    "allelically_excluded_gene",
    "gene_rearranged_at_DNA_level",
    "maternally_imprinted_gene",
    "paternally_imprinted_gene",
    "foreign_gene",
    "fusion_gene",
    "gene_cassette",
    "gene_with_non_canonical_start_codon",
    "gene_with_start_codon_CUG",
    "gene_with_polycistronic_transcript",
    "gene_with_dicistronic_transcript",
    "gene_with_dicistronic_mRNA",
    "gene_with_dicistronic_primary_transcript",
    "gene_with_trans_spliced_transcript",
    "mt_gene",
    "kinetoplast_gene",
    "maxicircle_gene",
    "minicircle_gene",
    "ncRNA_gene",
    "gRNA_gene",
    "lincRNA_gene",
    "miRNA_gene",
    "piRNA_gene",
    "RNase_MRP_RNA_gene",
    "RNase_P_RNA_gene",
    "rRNA_gene",
    "scRNA_gene",
    "snoRNA_gene",
    "snRNA_gene",
    "SRP_RNA_gene",
    "telomerase_RNA_gene",
    "tmRNA_gene",
    "tRNA_gene",
    "negatively_autoregulated_gene",
    "nuclear_gene",
    "nucleomorph_gene",
    "plasmid_gene",
    "plastid_gene",
    "apicoplast_gene",
    "chromoplast_gene",
    "ct_gene",
    "cyanelle_gene",
    "leucoplast_gene",
    "proplastid_gene",
    "positively_autoregulated_gene",
    "post_translationally_regulated_gene",
    "predicted_gene",
    "protein_coding_gene",
    "gene_with_edited_transcript",
    "gene_with_mRNA_with_frameshift",
    "gene_with_polyadenylated_mRNA",
    "gene_with_recoded_mRNA",
    "gene_with_mRNA_recoded_by_translational_bypass",
    "gene_with_stop_codon_read_through",
    "gene_with_stop_codon_redefined_as_pyrrollysine",
    "gene_with_stop_codon_redefined_as_selenocysteine",
    "gene_with_transcript_with_translational_frameshift",
    "proviral_gene",
    "endogenous_retroviral_gene",
    "non_functional_homolog_of_pseudogene",
    "non_processed_pseudogene",
    "cassette_pseudogene",
    "duplicated_pseudogene",
    "nuclear_mt_pseudogene",
    "pseudogene_by_unequal_crossing_over",
    "unitary_pseudogene",
    "polymorphic_pseudogene",
    "processed_pseudogene",
    "transposable_element_pseudogene",
    "recombinatorially_rearranged_gene",
    "recombinatorially_inverted_gene",
    "recombinatorially_rearranged_vertebrate_immune_system_gene",
    "rescue_gene",
    "wild_type_rescue_gene",
    "retrogene",
    "silenced_gene",
    "gene_silenced_by_DNA_modification",
    "gene_silenced_by_DNA_methylation",
    "gene_silenced_by_histone_modification",
    "gene_silenced_by_histone_deacetylation",
    "gene_silenced_by_histone_methylation",
    "gene_silenced_by_RNA_interference",
    "transgene",
    "floxed_gene",
    "translationally_regulated_gene",
    "transposable_element_gene",
    "engineered_foreign_transposable_element_gene"
}
"""GFF3 gene types as annotated by `SO 2.5.3 <http://www.sequenceontology.org/resources/intro.html>`_"""

_DEFAULT_GFF3_TRANSCRIPT_TYPES = {
     "transcript",
    "mature_transcript",
    "enzymatic_RNA",
    "mRNA",
    "mRNA_with_frameshift",
    "mRNA_with_minus_1_frameshift",
    "mRNA_with_minus_2_frameshift",
    "mRNA_with_plus_1_frameshift",
    "mRNA_with_plus_2_frameshift",
    "polyadenylated_mRNA",
    "polycistronic_mRNA",
    "dicistronic_mRNA",
    "monocistronic_mRNA",
    "recoded_mRNA",
    "mRNA_recoded_by_codon_redefinition",
    "mRNA_recoded_by_translational_bypass",
    "trans_spliced_mRNA",
    "ncRNA",
    "antisense_RNA",
    "MicF_RNA",
    "class_I_RNA",
    "class_II_RNA",
    "enhancerRNA",
    "guide_RNA",
    "lnc_RNA",
    "antisense_lncRNA",
    "intronic_lncRNA",
    "lincRNA",
    "piRNA",
    "priRNA",
    "rasiRNA",
    "RNase_MRP_RNA",
    "RNase_P_RNA",
    "rRNA",
    "large_subunit_rRNA",
    "rRNA_21S",
    "rRNA_23S",
    "rRNA_25S",
    "rRNA_28S",
    "rRNA_5_8S",
    "rRNA_5S",
    "small_subunit_rRNA",
    "rRNA_16S",
    "rRNA_18S",
    "rRNA_cleavage_RNA",
    "scRNA",
    "shRNA",
    "siRNA",
    "small_regulatory_ncRNA",
    "CsrB_RsmB_RNA",
    "DsrA_RNA",
    "GcvB_RNA",
    "IoR",
    "miRNA",
    "moR",
    "OxyS_RNA",
    "RNA_6S",
    "RprA_RNA",
    "RRE_RNA",
    "spot42_RNA",
    "tmRNA",
    "snoRNA",
    "C_D_box_snoRNA",
    "methylation_guide_snoRNA",
    "U14_snoRNA",
    "U3_snoRNA",
    "H_ACA_box_snoRNA",
    "pseudouridylation_guide_snoRNA",
    "snRNA",
    "U11_snRNA",
    "U12_snRNA",
    "U1_snRNA",
    "U2_snRNA",
    "U4_snRNA",
    "U4atac_snRNA",
    "U5_snRNA",
    "U6_snRNA",
    "U6atac_snRNA",
    "SRP_RNA",
    "tasiRNA",
    "telomerase_RNA",
    "telomeric_transcript",
    "anti_ARRET",
    "ARIA",
    "ARRET",
    "TERRA",
    "tRNA",
    "alanyl_tRNA",
    "arginyl_tRNA",
    "asparaginyl_tRNA",
    "aspartyl_tRNA",
    "cysteinyl_tRNA",
    "glutaminyl_tRNA",
    "glutamyl_tRNA",
    "glycyl_tRNA",
    "histidyl_tRNA",
    "isoleucyl_tRNA",
    "leucyl_tRNA",
    "lysyl_tRNA",
    "methionyl_tRNA",
    "phenylalanyl_tRNA",
    "prolyl_tRNA",
    "pyrrollysyl_tRNA",
    "selenocysteinyl_tRNA",
    "seryl_tRNA",
    "threonyl_tRNA",
    "tryptophanyl_tRNA",
    "tyrosyl_tRNA",
    "valyl_tRNA",
    "vault_RNA",
    "Y_RNA",
    "edited_transcript",
    "edited_mRNA",
    "edited_transcript_by_A_to_I_substitution",
    "non_functional_homolog_of_pseudogenic_transcript",
    "pseudogenic_transcript",
 }
"""GFF3 mature transcript types as annotated by `SO 2.5.3 <http://www.sequenceontology.org/resources/intro.html>`_"""

_DEFAULT_GFF3_EXON_TYPES = {
    "exon",
    "coding_exon",
    "noncoding_exon",
    "exon_of_single_exon_gene",
    "interior_exon",
    "interior_coding_exon",
    "five_prime_coding_exon",
    "three_prime_coding_exon"
    "five_prime_noncoding_exon",
    "three_prime_noncoding_exon",
    "pseudogenic_exon",
}
"""GFF3 exon feature types as annotated by `SO 2.5.3 <http://www.sequenceontology.org/resources/intro.html>`_"""

_DEFAULT_GFF3_CDS_TYPES = { 
     "CDS",
     "CDS_fragment",
     "CDS_indpendently_known",
     "CDS_predicted",
}
"""GFF3 CDS feature types as annotated by `SO 2.5.3 <http://www.sequenceontology.org/resources/intro.html>`_"""


#===============================================================================
# INDEX: Readers for GFF formats
#===============================================================================

StopFeature = SegmentChain(GenomicSegment("Stop",0,1,"."),type="StopFeature",ID="StopFeature")
"""Special |SegmentChain| emitted from GFF readers when the special line "###" is 
encountered, indicating that all previously returned features may be assembled 
into full objects. Also emitted when a GFF is sorted by chromosome, and
the chromosome name changes"""

class AbstractGFF_Reader(AbstractReader):
    """Abstract base class for GFF readers.
    
    Parses GFF streams line by line into |GenomicSegment|
    
    Attributes
    ----------
    metadata : dict
        Dictionary of metadata found in file headers
    """
    
    def __init__(self,stream,adjust_to_0=True,end_included=True,return_stopfeatures=True,is_sorted=False,tabix=False):
        """Create an |AbstractGFF_Reader|
        
        Parameters
        ----------
        stream : file-like
            Input stream pointing to GFF information
        
        adjust_to_0 : bool, optional
            Boolean, whether or not to adjust feature
            indices to a 0 base. True for GTF2 and GFF3
            files, as these are 1-indexed. (Default: True)
                                
        end_included : bool, optional
            Boolean, whether the end coordinate is
            included in the feature (closed or 'end-included' intervals)
            or not (half-open intervals). (Default: True)

        return_stopfeatures : bool, optional
            If True, will return a special |SegmentChain| called :py:obj:`StopFeature`
            signifying that all previously emitted GFF entries may be assembled
            into complete entities. These are emitted when the line "###"
            is encountered in a GFF. (Default: True)

        is_sorted : bool, optional
            If True and ``return_stopfeatures`` is True, assume the GFF is sorted.
            The reader will return :py:obj:`StopFeature` when the chromosome name
            of a given feature differs from that of the previous feature.
            (Default: False)
            
        tabix : boolean, optional
            Set to *True* if incoming streams are tabix-compressed, and
            using the parser :py:class:`pysam.asTuple` (Default: False)            
        """
        if tabix == True:
            stream = ("\t".join(X) for X in stream)
        
        self.chromosomes  = {}
        self.metadata   = {}
        self.adjust_to_0  = adjust_to_0
        self.end_included = end_included
        self.line_queue   = []
        self.return_stopfeatures = return_stopfeatures
        self.is_sorted = is_sorted
        
        line = next(stream)
        while line[0:2] == "##":
            self._parse_metatokens(line[2:])
            line = next(stream)
        
        self.line_queue.append(line)
        self._last_chrom = None
        
        self.stream = itertools.chain(self.line_queue,SkipBlankReader(stream))
        super(AbstractGFF_Reader,self).__init__(self.stream)

    def _parse_metatokens(self,inp):
        """Parses metadata embedded in a GFF stream, and stores
        these in appropriate attributes.
        
        Parameters
        ----------
        inp : str
            line of GFF input
            
        Raises
        ------
        StopIteration : when no features remain in file
        """
        items = inp.rstrip().split()
        if len(items) > 0:
            key = items[0]
            if key == "FASTA":
                raise StopIteration() #e.g. is end of features
            elif key == "sequence-region":
                try:
                    self.chromosomes[items[1]] = (items[2],items[3])
                except IndexError:
                    self.chromosomes[items[1]] = tuple(items[1:])
            elif key in self.metadata.keys():
                self.metadata[key] += ";" + " ".join(items[1:])
            else:
                self.metadata[key] = " ".join(items[1:])
    
    @abstractmethod
    def _parse_tokens(self,attr_string):
        """Placeholder function to parse column 9, which is formatted 
        differently in different GFF subtypes. Implement this 
        in subclasses

        Parameters
        ----------
        attr_string : str
            Ninth column of GFF

        Returns
        -------
        dict
            Dictionary of parsed tokens from ninth GFF column
        """
        pass
    
    def _parse_genomic_feature(self,line):
        """Parse GFF lines into |SegmentChain| objects

        Parameters
        ----------
        line : str
            Valid line of a GFF formatted file

        Returns
        -------
        |SegmentChain|
        """
        items = line.rstrip("\n").split("\t")
        chrom        = items[0]
        source       = items[1]
        feature_type = items[2]
        start        = int(items[3]) - int(self.adjust_to_0)
        end          = int(items[4]) - int(self.adjust_to_0) + int(self.end_included)
        score        = items[5]
        strand       = items[6]
        phase        = items[7]
        attr_string  = items[8]            
        info_dict = self._parse_tokens(attr_string)
        info_dict['source'] = source
        info_dict['score']  = score
        info_dict['phase']  = phase
        my_iv = GenomicSegment(chrom,start,end,strand)            
        info_dict['type'] = feature_type        
        my_feature = SegmentChain(my_iv,**info_dict)

        if chrom != self._last_chrom:
            old_chrom = self._last_chrom
            self._last_chrom = chrom
            if old_chrom is not None:
                if self.is_sorted == True and self.return_stopfeatures == True:
                    self.stream = itertools.chain([line],self.stream)
                    return StopFeature

        return my_feature    
    
    def filter(self,line):
        """Parses lines of the GFF3 stream into |SegmentChain|
        When metadata is found, temporarily delegates processing to 
        :meth:`_parse_metatokens`, and then reads the next genomic feature
        
        Parameters
        ----------
        line
            Next line from GFF3 stream
        
        Returns
        -------
        |SegmentChain|
            Next feature in file
        """
        if line[0:3] == "###":
            if self.return_stopfeatures == True:
                return StopFeature
            else:
                return self.__next__()
        elif line[0:2] == "##":
            self._parse_metatokens(line[2:])
            return self.__next__()
        elif line[0:1] == "#":
            return self.__next__()
        else:
            return self._parse_genomic_feature(line)

    
class GFF3_Reader(AbstractGFF_Reader):
    """Parses `GFF3`_ streams line by line into |SegmentChain|
    Assumes input stream to be GFF3-compliant, as specified at
    the `Sequence Ontology GFF3 specification <http://song.sourceforge.net/gff3.shtml>`_.
    
    In short, this means the file is 1-indexed, uses a controlled
    vocabulary, and follows a defined schema of parent-child
    relationships between features. It is unclear whether
    feature coordinates are end-included or not, though several 
    organizations (Flybase, SGD, IGV) use end-included GFFs.
    
    1-based coordinates are adjusted to a 0 base in keeping with Python
    conventions.
    
    GFF3 attributes (from column 9) are stored in a dictionary called ``attr``
    in each of the returned |SegmentChain| objects as unescaped strings. The values
    for the attributes 'Parent', 'Alias', 'Dbxref', 'dbxref', and 'Note',
    if present, are lists rather than strings, because the `GFF3`_ spec enables
    these to have multiple values. 
    
    Attributes
    ----------
    metadata : dict
        Dictionary of metadata found in file headers
    """
    
    def __init__(self,stream,end_included=True,return_stopfeatures=False,is_sorted=False,tabix=False):
        """Create a |GFF3_Reader|
        
        Parameters
        ----------
        stream : file-like
            Input stream pointing to GFF3 information
        
        end_included : bool, optional
            Boolean, whether the end coordinate is
            included in the feature (closed or 'end-included' intervals)
            or not (half-open intervals). (Default: True)

        return_stopfeatures : bool, optional
            If True, will return a special |SegmentChain| called :py:obj:`StopFeature`
            signifying that all previously emitted GFF entries may be assembled
            into complete entities. These are emitted when the line "###"
            is encountered in a GFF. (Default: False)

        is_sorted : bool, optional
            If True and ``return_stopfeatures`` is True, assume the GFF is sorted.
            The reader will return :py:obj:`StopFeature` when the chromosome name
            of a given feature differs from that of the previous feature.
            (Default: False)
            
        tabix : boolean, optional
            Set to *True* if incoming streams are tabix-compressed, and
            using the parser :py:class:`pysam.asTuple` (Default: False)            
         """
        super(GFF3_Reader,self).__init__(stream,
                                         adjust_to_0=True,
                                         end_included=end_included,
                                         return_stopfeatures=return_stopfeatures,
                                         is_sorted=is_sorted,tabix=tabix)
    
    def _parse_tokens(self,inp):
        """Parse column 9 of GFF3 into dictionary

        Parameters
        ----------
        inp : str
            Ninth column of GFF

        Returns
        -------
        dict
            Dictionary of parsed tokens from ninth GFF3 column
        """
        return parse_GFF3_tokens(inp)
    
            
class GTF2_Reader(AbstractGFF_Reader):
    """Parses GTF2 streams line by line into |SegmentChain|
    Assumes input stream to be GTF2-compliant, as specified at
    the `GTF2 file specification <http://mblab.wustl.edu/GTF22.html>`_
    
    In short, this means the file is 1-indexed, and every feature has
    "gene_id" and "transcript_id" attributes.
    
    1-based coordinates are adjusted to a 0 base in keeping with Python
    conventions.
    
    Attributes
    ----------
    metadata : dict
        Dictionary of metadata found in file headers    
    """
    def __init__(self,stream,end_included=True,return_stopfeatures=False,is_sorted=False,tabix=False):
        """Create a |GTF2_Reader|
        
        Parameters
        ----------
        stream : file-like
            Input stream pointing to GFF3 information
        
        end_included : bool, optional
            Boolean, whether the end coordinate is
            included in the feature (closed or 'end-included' intervals)
            or not (half-open intervals). (Default: True)

        return_stopfeatures : bool, optional
            If True, will return a special |SegmentChain|  called :py:obj:`StopFeature`
            signifying that all previously emitted GFF entries may be assembled
            into complete entities. These are emitted when the line "###"
            is encountered in a GTF2. (Default: False)

        is_sorted : bool, optional
            If True and ``return_stopfeatures`` is True, assume the GTF2 is sorted
            by chromosome. The reader will return :py:obj:`StopFeature` when
            the chromosome name of a given feature differs from that of the previous
            feature.
            (Default: False)
            
        tabix : boolean, optional
            Set to *True* if incoming streams are tabix-compressed, and
            using the parser :py:class:`pysam.asTuple` (Default: False)            
         """
        super(GTF2_Reader,self).__init__(stream,
                                         adjust_to_0=True,
                                         end_included=end_included,
                                         return_stopfeatures=return_stopfeatures,
                                         is_sorted=is_sorted,tabix=tabix)
        
    def _parse_tokens(self,inp):
        """Parse column 9 of GTF2 into dictionary

        Parameters
        ----------
        inp : str
            Ninth column of GFF

        Returns
        -------
        dict
            Dictionary of parsed tokens from ninth GTF2 column
        """
        return parse_GTF2_tokens(inp)


class AbstractGFF_Assembler(AssembledFeatureReader):
    """Abstract base class for readers that assemble composite features from one
    or more features in one or more streams of GTF2 or GFF3 data.
    
    
    Attributes
    ----------
    stream : file-like
        Input stream, usually constructed from or more open filehandles
    
    metadata : dict
        Various attributes gleaned from the stream, if any

    counter : int
        Cumulative line number counter over all streams

    printer : file-like, optional
        Logger implementing a ``write()`` method.

    rejected : list
        A list of transcript IDs that failed to assemble properly
    """
    
    def __init__(self,*streams,**kwargs):
        """Create a |AbstractGFF_Assembler|
        
        Parameters
        ----------
        streams : file-like
            One or more open streams of input data.

        is_sorted : bool, optional
            GFF is sorted by chromosome name, allowing some memory savings
            (Default: False)
        
        return_type : |SegmentChain| or subclass, optional
            Type of feature to return from assembled subfeatures (Default: |SegmentChain|)

        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations. If set to
            True, three nucleotides will be added to the threeprime end of each
            CDS annotation. Default: False
        
        printer : file-like, optional
            Logger implementing a ``write()`` method. Default: |NullWriter|
        
        reader_class : class
            GFF3_Reader or GTF2_Reader

        tabix : boolean, optional
            Set to *True* if incoming streams are tabix-compressed, and
            using the parser :py:class:`pysam.asTuple` (Default: False)
                    
        **kwargs
            Other keyword arguments used by specific parsers
        """
        
        tabix = kwargs.get("tabix",False)
        end_included = kwargs.get("end_included",True)
        return_stopfeatures = kwargs.get("return_stopfeatures",True)
        is_sorted = kwargs.get("is_sorted",False)
        reader_class = kwargs.get("reader_class")
        
        iterables = []
        for stream in streams:
            iterables.append(reader_class(stream,
                                          end_included=end_included,
                                          return_stopfeatures=return_stopfeatures,
                                          is_sorted=is_sorted,
                                          tabix=tabix))
            iterables.append([StopFeature])
            
        self.stream = itertools.chain.from_iterable(iterables)
        

        self.printer       = kwargs.get("printer",NullWriter())
        self.return_type   = kwargs.get("return_type",SegmentChain)        
        self.add_three_for_stop = kwargs.get("add_three_for_stop",False)

        self.metadata = {}
        self.rejected = []
        self.counter  = 0
        
        self._transcript_cache = iter([])
        self._feature_cache = {}
    
    def _finalize(self,tx):
        return tx
    
    def __iter__(self):
        return self
    
    def next(self):
        return self.__next__()

    @abstractmethod
    def _collect(self,feature):
        """Collect relevant features of transcripts
        
        Parameters
        ----------
        feature : |SegmentChain|
            Feature to collect
        """        
        pass
    
    @abstractmethod
    def _assemble_transcripts(self):
        """Assemble |Transcript| s
        
        Returns
        -------
        list
            list of transcripts
        """        
        pass
    
    @abstractmethod
    def _reset(self):
        """Release memory and reset internal hashes"""
        pass
        
    def __next__(self):
        """Return next transcript in GTF2/GFF3, using lazy evaluation as follows:
         
        #.  If there exist assembled transcripts in ``self._transript_cache``, 
            return the next transcript. Transcripts in the cache are stored
            lexically.
         
        #.  Otherwise, collect features from the GTF2/GFF3 stream until either
            a StopFeature or EOF is encountered. At that point, assemble transcripts
            and store them in the iterator ``self._transcript_cache``
         
        Returns
        -------
        |SegmentChain|
            Next complex feature in annotation (usually a transcript)
        """
        try:
            return self._finalize(next(self._transcript_cache))
        except StopIteration:
            # read GTF2/GFF3, populating feature_cache and building transcripts
            # only after termination triggers (StopFeature or EOF)
            try: # if GTF2/GFF3 features, collect until StopIteration
                for feature in self.stream:
                    self.counter += 1
                    if feature == StopFeature: # if GTF2 is sorted or contains ###, raise stop, assemble, and free memory
                        raise StopIteration()
                    else:
                        self._collect(feature)
            except StopIteration: # end of GTF2/GFF3 feature stream, assemble transcripts
                self.printer.write("Assembling next batch of transcripts...")
                transcripts, rejected = self._assemble_transcripts()
                if len(transcripts) > 0:
                    transcripts = sorted(transcripts,key=sort_segmentchains_lexically)
                    self._transcript_cache = iter(transcripts)
                    self.rejected.extend(rejected)
                    self._reset()
                else: # happens when we get two StopFeatures in a row
                    return self.__next__()
            return self._finalize(next(self._transcript_cache))


class GTF2_TranscriptAssembler(AbstractGFF_Assembler):
    """Assemble features in one or more streams of GTF2 data into an iterator over |Transcript| objects,
    collecting exon and CDS features based upon shared ``transcript_id``s.
    Attributes that have common values for all exons and CDS within a transcript
    are propagated up to the ``attr` dict of the assembled |Transcript|. Other
    attributes from individual CDS or exon components are discarded.
     
    Transcripts are returned in lexical order.


    Attributes
    ----------
    stream : file-like
        Input stream, usually constructed from or more open filehandles
    
    metadata : dict
        Various attributes gleaned from the stream, if any

    counter : int
        Cumulative line number counter over all streams

    printer : file-like, optional
        Logger implementing a ``write()`` method.

    rejected : list
        A list of transcript IDs from transcripts that failed to assemble properly


    Notes
    -----
    To save memory, transcripts are assembled using lazy evaluation.
    Assembly proceeds as follows:
    
        #.  If there exist assembled transcripts in ``self._transript_cache``, 
            return the next transcript. Transcripts in the cache are stored
            lexically.
        
        #.  Otherwise, collect features from the GTF2 stream until either a
            '###' line or EOF is encountered. Then, assemble transcripts and
            store them in ``self._transcript_cache``. Delete unused features
            from memory. If the GTF2 is sorted, then a change in chromosome
            name will also trigger assembly of collected features.
    """

    # transcripts can be represented as collections of exons + cds
    # or cds + UTRs, et c. We consider all UTR and exons as exons
    # and CDS, and start & stop codons as CDS areas
    _feature_map = { "exon"  : ["exon_like"],
                     "5UTR"  : ["exon_like"],
                     "3UTR"  : ["exon_like"],
                     "CDS"   : ["CDS_like","exon_like"],
                     "start_codon" : ["CDS_like","exon_like"],
                     "stop_codon"  : ["CDS_like","exon_like"],
                   }
    dtmp = { "exon_like" : {}, "CDS_like" : {} }

    def __init__(self,*streams,**kwargs):
        """Create a |GTF2_TranscriptAssembler|
        
        Parameters
        ----------
        streams : file-like
            One or more open streams of GTF2 data.

        is_sorted : bool, optional
            GTF2 is sorted by chromosome name, allowing some memory savings
            (Default: False)
        
        return_type : |SegmentChain| or subclass, optional
            Type of feature to return from assembled subfeatures (Default: |SegmentChain|)

        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations. If set to
            *True*, three nucleotides will be added to the threeprime end of each
            CDS annotation, UNLESS the annotated transcript contains explicit stop_codon 
            feature. (Default: False)
        
        printer : file-like, optional
            Logger implementing a ``write()`` method. Default: |NullWriter|
            
        tabix : boolean, optional
            Set to *True* if incoming streams are tabix-compressed, and
            using the parser :py:class:`pysam.asTuple` (Default: False)            
        """
        AbstractGFF_Assembler.__init__(self,*streams,reader_class=GTF2_Reader,**kwargs)
        self._feature_cache = { "exon_like" : {}, "CDS_like" : {}}


    def _collect(self,feature):
        """Collect relevant features of transcripts, and populate ``self._feature_cache`
        
        Parameters
        ----------
        feature : |SegmentChain|
            Feature to collect
        """
        feature_classes = self._feature_map.get(feature.attr["type"],None)
        if feature_classes is not None:
            tname = feature.attr.get("transcript_id")
            for feature_class in feature_classes:
                try:
                    self._feature_cache[feature_class][tname].append(feature)
                except KeyError:
                    self._feature_cache[feature_class][tname] = [feature]
                                    
    def _assemble_transcripts(self):
        """Assemble |Transcript| s from features in ``self._feature_cache``
        
        Returns
        -------
        list
            list of transcripts
        """
        transcripts, rejected = _assemble_transcripts_from_gtf2_dicts(self._feature_cache,
                                                                      printer=self.printer,
                                                                      add_three_for_stop=self.add_three_for_stop)
        return transcripts, rejected
        
    def _reset(self):
        """Release memory and reset internal hashes"""
        del self._feature_cache
        gc.collect()
        del gc.garbage[:]
        self._feature_cache = { "exon_like" : {}, "CDS_like" : {}}


class GFF3_TranscriptAssembler(AbstractGFF_Assembler):
    """Assemble features in one or more streams of GFF3 data into an iterator over |Transcript| objects,
    
    Transcripts are returned in lexical order.

    Attributes
    ----------
    stream : file-like
        Input stream, usually constructed from or more open filehandles
    
    metadata : dict
        Various attributes gleaned from the stream, if any

    counter : int
        Cumulative line number counter over all streams

    printer : file-like, optional
        Logger implementing a ``write()`` method.

    rejected : list
        A list of transcript IDs from transcripts that failed to assemble properly


    Notes
    -----
    GFF3 schemas vary
        Different GFF3s have different schemas of hierarchy. We deal with that here
        by allowing users to supply `transcript_types` and `exon_types`, to indicate
        which sorts of features should be included.

    Identity relationships between elements vary between GFF3 files
        Also, different GFF3s specify discontiguous features differently. For example,
        in Flybase, different exons of a transcript will have unique IDs, but will share
        the same "Parent" attribute in column 9 of the GFF. In Wormbase, however, different
        exons of the same transcript will share the same ID. Here, we treat GFFs as if
        they are written in the Flybase style. We may support alternate formats in the future.    
    
    Transcript assembly
        To save memory, transcripts are assembled using lazy evaluation.
        Assembly proceeds as follows:
    
        #.  If there exist assembled transcripts in ``self._transript_cache``, 
            return the next transcript. Transcripts in the cache are stored
            lexically.
        
        #.  Otherwise, collect features from the GFF3 stream until either a
            '###' line or EOF is encountered. Then, assemble transcripts and
            store them in ``self._transcript_cache``. Delete unused features
            from memory. If the GFF3 is sorted, then a change in chromosome
            name will also trigger assembly of collected features.
    """

    def __init__(self,*streams,**kwargs):
        """Create a |GFF3_TranscriptAssembler|
        
        Parameters
        ----------
        streams : file-like
            One or more open streams of GFF3 data.
        
        is_sorted : bool, optional
            GFF3 is sorted by chromosome name, allowing some memory savings
            (Default: False)
        
        return_type : |SegmentChain| or subclass, optional
            Type of feature to return from assembled subfeatures (Default: |SegmentChain|)

        add_three_for_stop : bool, optional
            Some annotation files exclude the stop codon from CDS annotations. If set to
            True, three nucleotides will be added to the threeprime end of each
            CDS annotation. Default: False
        
        transcript_types : list, optional
            List of GFF3 feature types that should be considered as transcripts
            (Default: as specified in SO 2.5.3 )

        exon_types : list, optional
            List of GFF3 feature types that should be considered as exons or
            contributing to transcript nucleotide positions
            during transcript assembly (Default: as specified in SO 2.5.3 )
        
        cds_types : list, optional
            List of GFF3 feature types that should be considered as CDS or
            contributing to transcript coding regions during transcript assembly
            (Default: as specified in SO 2.5.3 )
        
        printer : file-like, optional
            Logger implementing a ``write()`` method. Default: |NullWriter|
            
        tabix : boolean, optional
            Set to *True* if incoming streams are tabix-compressed, and
            using the parser :py:class:`pysam.asTuple` (Default: False)
        
        
        Notes
        -----
        Sequence Ontology 2.5.3
            By default, this assembler constructs transcripts following the GFF3
            schema from the `SO Consortium <http://www.sequenceontology.org/resources/intro.html>`_
        """
        AbstractGFF_Assembler.__init__(self,*streams,reader_class=GFF3_Reader,**kwargs)
        self.transcript_types = set(kwargs.get("transcript_types",_DEFAULT_GFF3_TRANSCRIPT_TYPES))
        self.exon_types       = set(kwargs.get("exon_types",_DEFAULT_GFF3_EXON_TYPES))
        self.cds_types = set(kwargs.get("cds_types",_DEFAULT_GFF3_CDS_TYPES))
        self.transcript_components = self.exon_types | self.cds_types
        self._feature_cache = {}
        self._tx_features = {}
        self._reset()

    def _collect(self,feature):
        """Collect relevant features of transcripts, and populate ``self._feature_cache`
        
        Parameters
        ----------
        feature : |SegmentChain|
            Feature to collect
        """
        feature_name = feature.get_name()
        
        if feature.attr["type"] in self.transcript_types:
            try:
                self._tx_features[feature_name].append(feature)
            except KeyError:
                self._tx_features[feature_name] = [feature]
        
        elif feature.attr["type"] in self.transcript_components:
            tnames = feature.attr.get("Parent")
            for tname in tnames:
                try:
                    self._feature_cache[feature.attr["type"]][tname].append(feature)
                except KeyError:
                    self._feature_cache[feature.attr["type"]][tname] = [feature]
                                    
    def _assemble_transcripts(self):
        """Assemble |Transcript| s from features in ``self._feature_cache``
        
        Returns
        -------
        list
            list of transcripts
        """
        rejected   = []
        transcripts = []
        tx_features_counted = []
        
        for type_ in self.exon_types:
            for tname in self._feature_cache[type_].keys():
                tx_features_counted.append(tname)
                exons = self._feature_cache[type_].get(tname,[])
                cds = []
                for cds_type in self.cds_types:
                    cds.extend(self._feature_cache[cds_type].get(tname,[]))
                
                gene_id = self._tx_features[tname][0].attr.get("Parent",tname) # use transcript name as gene if no Parent
                # gene IDs are now returned as lists from GFF3 parser
                gene_id = ",".join(sorted(gene_id))

                attr    = self._tx_features[tname][0].attr #get attr from transcript object
                attr["ID"] = tname
                attr["transcript_id"] = tname
                attr["gene_id"] = gene_id
                exon_segments = [X.spanning_segment for X in exons]                        

                if len(cds) > 0:
                    cds   = sorted(cds,key = lambda x: x.spanning_segment.start)
                    attr["cds_genome_start"] = cds[0].spanning_segment.start
                    attr["cds_genome_end"] = cds[-1].spanning_segment.end
        
                    exons = sorted(exons,key = lambda x: x.spanning_segment.start)
                    # correct exon boundaries that don't include entire CDS
                    if cds[0].spanning_segment.start < exons[0].spanning_segment.start:
                        exons[0].spanning_segment.start = cds[0].spanning_segment.start
                    if cds[-1].spanning_segment.end > exons[-1].spanning_segment.end:
                        exons[-1].spanning_segment.end = cds[-1].spanning_segment.end
                try:
                    my_tx = Transcript(*tuple(exon_segments),**attr)
                    if self.add_three_for_stop == True:
                        my_tx = add_three_for_stop_codon(my_tx)
                    
                    transcripts.append(my_tx)     
                                   
                except AssertionError:
                    self.printer.write("Rejecting %s because it contains exons on multiple strands." % tname)
                    # transcripts with exons on two strands
                    rejected.append(tname)
                except KeyError:
                    # transcripts where CDS ends outside bounds of transcript
                    # there are 25 of these in flybase r5.43
                    self.printer.write("Rejecting %s because start or stop codons are outside exon boundaries." % tname)                        
                    rejected.append(tname)
            
        tx_features_not_counted = set(self._tx_features.keys()) - set(tx_features_counted)
        for txid in tx_features_not_counted:
            attr = self._tx_features[txid][0].attr
            attr["ID"] = txid
            attr["transcript_id"] = txid
            segments  = [X.spanning_segment for X in self._tx_features[txid]]
            my_txmodel = Transcript(*tuple(segments),**attr)
            transcripts.append(my_txmodel)                    
        
        return transcripts, rejected
        
    def _reset(self):
        """Release memory and reset internal hashes"""
        del self._feature_cache
        del self._tx_features
        gc.collect()
        del gc.garbage[:]
        self._tx_features = {}  
        self._feature_cache = { X : copy.deepcopy({}) for X in self.transcript_components }
            

        

#===============================================================================
# INDEX: convenience functions
#===============================================================================
@deprecated
def GFF3_to_Transcripts(stream,end_included=True,
                        transcript_types=_DEFAULT_GFF3_TRANSCRIPT_TYPES,
                        exon_types=_DEFAULT_GFF3_EXON_TYPES,
                        add_three_for_stop=False,
                        printer=None):
    """Parses a GFF3 stream into a list of assembled |Transcript|
    
    Parameters
    ----------
    stream : file-like
        Input stream pointing to GFF3 information
    
    end_included : bool, optional
        Boolean, whether the end coordinate is
        included in the feature (closed or 'end-included' intervals)
        or not (half-open intervals). (Default: *True*)
    
    printer : file-like, optional
        Filehandle or sys.stderr-like for logging
            
    transcript_types : list<str>, optional
        Types of feature to include as transcripts, even if no exons
        or CDS.

    exon_types : list<str>, optional
        Types of feature to include as exons
        (e.g. *["exon","noncoding_exon"]*, et c; Default: *["exon"]*)
    

    add_three_for_stop : bool, optional
        Some GFF3 files exclude the stop codon from CDS annotations. If set to
        True, three nucleotides will be added to the threeprime end of each
        CDS annotation. Default: *False*
        
    Returns
    -------
    list<|Transcript|>
        Assembled transcripts
        
    list<transcript>
        Transcripts that were malformed in annotaiton file or otherwise rejected by parser
    
    
    Notes
    -----
    GFF3 schemas vary
        Different GFF3s have different schemas of hierarchy. We deal with that here
        by allowing users to supply `transcript_types` and `exon_types`, to indicate
        which sorts of features should be included.

    Identity relationships vary
        Also, different GFF3s specify discontiguous features differently. For example,
        in Flybase, different exons of a transcript will have unique IDs, but will share
        the same "Parent" attribute in column 9 of the GFF. In Wormbase, however, different
        exons of the same transcript will share the same ID. Here, we treat GFFs as if
        they are written in the Flybase style. We may support alternate formats in the future.
    """
    tx_components = set(exon_types) | set(["CDS"])
    tx_features = {}
    dtmp = { "exon" : {}, "CDS" : {} }
    transcripts_out = {}
    rejected_transcripts = []
    tx_features_counted  = []
    for feature in GFF3_Reader(stream,end_included=end_included):
        if feature.attr["type"] in tx_components: #dtmp.keys():
            
            if feature.attr["type"] in exon_types:
                feature.attr["type"] = "exon"
                
            tnames = feature.attr.get("Parent",[])
            for tname in tnames:
                try:
                    dtmp[feature.attr["type"]][tname].append(feature)
                except KeyError:
                    dtmp[feature.attr["type"]][tname] = [feature]
        elif feature.attr["type"] in transcript_types:
            try:
                tx_features[feature.get_name()].append(feature)
            except KeyError:
                tx_features[feature.get_name()] = [feature]
                
    for tname in dtmp["exon"].keys():
        tx_features_counted.append(tname)
        exons = dtmp["exon"].get(tname,[])
        cds = dtmp["CDS"].get(tname,[])
        gene_id = tx_features[tname][0].attr.get("Parent",tname) # use transcript name as gene if no Parent
        gene_id = ",".join(sorted(gene_id))
        attr    = tx_features[tname][0].attr #get attr from transcript object
        attr["ID"] = tname
        attr["transcript_id"] = tname
        attr["gene_id"] = gene_id
        exon_segments = [X.spanning_segment for X in exons]
        if len(cds) > 0:
            cds   = sorted(cds,key = lambda x: x.spanning_segment.start)
            attr["cds_genome_start"] = cds[0].spanning_segment.start
            attr["cds_genome_end"] = cds[-1].spanning_segment.end

            exons = sorted(exons,key = lambda x: x.spanning_segment.start)
            # correct exon boundaries that don't include entire CDS
            if cds[0].spanning_segment.start < exons[0].spanning_segment.start:
                exons[0].spanning_segment.start = cds[0].spanning_segment.start
            if cds[-1].spanning_segment.end > exons[-1].spanning_segment.end:
                exons[-1].spanning_segment.end = cds[-1].spanning_segment.end
        try:    
            my_tx = Transcript(*tuple(exon_segments),**attr)
            if add_three_for_stop == True:
                my_tx = add_three_for_stop_codon(my_tx)
            transcripts_out[tname] = my_tx
        except AssertionError:
            if printer is not None:
                printer.write("Rejecting %s because it contains exons on multiple strands." % tname)
            # transcripts with exons on two strands
            rejected_transcripts.append(tname)
        except KeyError:
            # transcripts where CDS ends outside bounds of transcript
            # there are 25 of these in flybase r5.43
            rejected_transcripts.append(tname)
            if printer is not None:
                printer.write("Rejecting %s because start or stop codons are outside exon boundaries." % tname)
    
    # which transcripts did not have exons or CDS? 
    # this occurs in GFFs where exons are not explicitly annotated (e.g. yeast pseudogenes, et c)
    # n.b. these will not be double-counted, because we're ignoring those we have already taken
    tx_features_not_counted = set(tx_features.keys()) - set(tx_features_counted)
    for txid in tx_features_not_counted:
        attr = tx_features[txid][0].attr
        attr["ID"] = txid
        attr["transcript_id"] = txid
        segments  = [X.spanning_segment for X in tx_features[txid]]
        my_txmodel = Transcript(*tuple(segments),**attr)
        transcripts_out[txid] = my_txmodel
            
    return sorted(transcripts_out.values(),key=sort_segmentchains_lexically), list(set(rejected_transcripts))

@deprecated
def GTF2_to_Transcripts(stream,end_included=True,printer=NullWriter(),add_three_for_stop=False,is_sorted=False):
    """Parses a GTF2 stream into a list of assembled |Transcript|
    
    Parameters
    ----------
    stream : file-like
        Input stream pointing to GTF2 information
    
    end_included : bool
        Boolean, whether the end coordinate is
        included in the feature (closed or 'end-included' intervals)
        or not (half-open intervals). (Default: True)
    
    printer : file-like
        Filehandle or sys.stderr-like for logging
        
    add_three_for_stop : bool, optional
        Some annotation files exclude the stop codon from CDS annotations. If set to
        *True*, three nucleotides will be added to the threeprime end of each
        CDS annotation, UNLESS the annotated transcript contains explicit stop_codon 
        feature. (Default: False)
    
    is_sorted : bool, optional
        If True, assume GTF2 is sorted by chromosome. In this case, transcripts
        will be assembled from component objects, and memory emptied every
        time the chromosome name changes.

    Returns
    -------
    list<|Transcript|>
        Assembled transcripts
        
    list<transcript>
        Transcripts that were malformed in annotaiton file or otherwise rejected by parser
    """
    # transcripts can be represented as collections of exons + cds
    # or cds + UTRs, et c. We consider all UTR and exons as exons
    # and CDS, and start & stop codons as CDS areas
    feature_map = { "exon"  : ["exon_like"],
                     "5UTR"  : ["exon_like"],
                     "3UTR"  : ["exon_like"],
                     "CDS"   : ["CDS_like","exon_like"],
                     "start_codon" : ["CDS_like","exon_like"],
                     "stop_codon"  : ["CDS_like","exon_like"],
                   }
    dtmp = { "exon_like" : {}, "CDS_like" : {} }
    transcripts = []
    rejected_transcripts = []
    for n,feature in enumerate(GTF2_Reader(stream,end_included=end_included,is_sorted=is_sorted)):
        if feature == StopFeature:
            printer.write("Assembling next batch of %s features..." % n)
            tr, rej = _assemble_transcripts_from_gtf2_dicts(dtmp,printer=printer,add_three_for_stop=add_three_for_stop)
            transcripts.extend(tr)
            rejected_transcripts.extend(rej)
            del dtmp
            gc.collect()
            del gc.garbage[:]
            dtmp = { "exon_like" : {}, "CDS_like" : {} }
        else:
            feature_classes = feature_map.get(feature.attr["type"],None)
            if feature_classes is not None:
                tname = feature.attr.get("transcript_id")
                for feature_class in feature_classes:
                    try:
                        dtmp[feature_class][tname].append(feature)
                    except KeyError:
                        dtmp[feature_class][tname] = [feature]

    tr, rej = _assemble_transcripts_from_gtf2_dicts(dtmp,printer=printer,add_three_for_stop=add_three_for_stop)
    transcripts.extend(tr)
    rejected_transcripts.extend(rej)
            
    return sorted(transcripts,key=sort_segmentchains_lexically), list(set(rejected_transcripts))

def _assemble_transcripts_from_gtf2_dicts(cds_exon,add_three_for_stop=False,printer=NullWriter()):
    """Assemble |Transcript| objects from a dictionary of features
    mapping transcript IDs to corresponding CDS and exon features.
    
    Attributes common to all CDS and exons for a given transcript (e.g. 
    ``gene_id`` and ``transcript_id``) are propagated up to the |Transcript|.
    Other component attributes are discarded.
    
    Parameters
    ----------
    cds_exon : dict
        Dictionary of dictionaries with keys "CDS_like" and "exon_like" that map
        transcript IDs to their CDS and exon components (from |GTF2_Reader|)
    
    add_three_for_stop : bool, optional
        Some annotation files exclude the stop codon from CDS annotations. If set to
        *True*, three nucleotides will be added to the threeprime end of each
        CDS annotation, UNLESS the annotated transcript contains explicit stop_codon 
        feature. (Default: False)
    
    printer : file-like, optional
        Logger implementing a ``write()`` method. Default: |NullWriter|
    """
    rejected_transcripts = []
    transcripts = []
    for tname in set(cds_exon["exon_like"].keys()) | set(cds_exon["CDS_like"].keys()):
        exons = cds_exon["exon_like"].get(tname,[])
        cds = cds_exon["CDS_like"].get(tname,[])
        if len(exons) > 0:
            #gene_id = exons[0].attr.get("gene_id")
            #attr    = exons[0].attr
            exons = sorted(exons,key = lambda x: x.spanning_segment.start)
            exon_segments = [X.spanning_segment for X in exons]
        elif len(cds) > 0:
            # if cds but no exons, create exons since they are implied
            #gene_id = cds[0].attr.get("gene_id")
            #attr    = cds[0].attr
            exons = sorted(cds,key = lambda x: x.spanning_segment.start)
            exon_segments = [X.spanning_segment for X in exons]
        
        # propagate attributes that are the same in all exons/cds
        # to parent. This should include `gene_id` and `transcript_id`
        attr = get_identical_attributes(exons + cds)
        #attr["ID"] = tname
        #attr["transcript_id"] = tname
        #attr["gene_id"] = gene_id        
        if len(cds) > 0:
            cds = sorted(cds,key = lambda x: x.spanning_segment.start)
            attr["cds_genome_end"]   = cds[-1].spanning_segment.end
            attr["cds_genome_start"] = cds[0].spanning_segment.start
        
        try:
            my_tx = Transcript(*tuple(exon_segments),**attr)
            if add_three_for_stop == True:
                # only move stop codons if no exon feature is of type "stop_codon"
                if "stop_codon" not in set([X.attr["type"] for X in exons]): 
                    my_tx = add_three_for_stop_codon(my_tx)
                    
            transcripts.append(my_tx)
        except AssertionError:
            printer.write("Rejecting %s because it contains exons on multiple strands." % tname)
            # transcripts with exons on two strands
            rejected_transcripts.append(tname)
        except KeyError:
            # transcripts where CDS ends outside bounds of transcript
            # there are 25 of these in flybase r5.43
            rejected_transcripts.append(tname)
            printer.write("Rejecting %s because start or stop codons are outside exon boundaries." % tname)

    return sorted(transcripts,key=sort_segmentchains_lexically), rejected_transcripts
