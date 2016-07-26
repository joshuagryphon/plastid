#!/usr/bin env python
"""Support for reading features from Genbank (.gb) files.
"""
import gc

from Bio import SeqIO
from plastid.genomics.roitools import (SegmentChain,
                                       Transcript,
                                       add_three_for_stop_codon
                                      )
from plastid.readers.gff import StopFeature 
from plastid.readers.gff import AbstractGFF_Assembler
from plastid.util.io.openers import multiopen, NullWriter
from plastid.util.services.decorators import notimplemented


# actually a function but it behaves like the rest of our reader classes
# internals of BioPython made it easier to code this way anyway
def GenbankReader(*streams,**kwargs):
    """Read individual features from Genbank files
    ]
    Parameters
    ----------
    *streams
    
    return_stopfeatures : bool, optional
        If `True`, will return a special |SegmentChain|  called :py:obj:`StopFeature`
        signifying that all previously emitted SegmentChains may be assembled
        into complete entities. These are emitted when the contigs change.
    """
    return_stopfeatures = kwargs.get("return_stopfeatures",False)
    files = multiopen(streams,fn=open)

    for file_ in files:
        gbstream = SeqIO.parse(file_,"genbank")
        for seq in gbstream:
            for feature in seq.features:
                yield SegmentChain.from_biopython_seqfeature(feature, seq.id)
                
            if return_stopfeatures == True:
                yield StopFeature



# See
#   - http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/
#   - http://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/#protein_id
# locus_tag -> gene_id
# gene -> Name
# for eukaryotes, CDS have transcript_id that match mRNA transcript_id 
# gene of CDS and mRNA should match gene or locus tag
# all CDS have product, 
# trnascript_id may be protein_id + 't' or 'mrna'
@notimplemented
class GenbankTranscriptAssembler(AbstractGFF_Assembler):
    """Assemble transcripts from features in Genbank files
    """
    _feature_map = { "exon"  : ["exon_like"],
                     "mRNA"  : ["exon_like"],
                     "CDS"   : ["CDS_like","exon_like"],
                     "start_codon" : ["CDS_like","exon_like"],
                     "stop_codon"  : ["CDS_like","exon_like"],
                   }
    
    def __init__(self,*streams,**kwargs):
        kwargs["return_stopfeature"] = True
        
        self.reader  = GenbankReader(*streams,**kwargs)                   
        self.counter = 0
        self.printer = kwargs.get("printer",NullWriter())

        self.return_type   = kwargs.get("return_type",SegmentChain)        
        add_three_for_stop = kwargs.get("add_three_for_stop",False)
        self._finalize =  add_three_for_stop_codon if add_three_for_stop == True else lambda x: x

        self.metadata    = {}
        self.rejected    = []
        self._transcript_cache = iter([])
        self._feature_cache = { "mrna_like" : {}, "CDS_like" : {}}

    def _collect(self,feature):
        """Collect relevant features of transcripts
        
        Parameters
        ----------
        feature : |SegmentChain|
            Feature to collect
        """        
        pass
    
    def _assemble_transcripts(self):
        """Assemble |Transcript| objects from collected features 
        
        Returns
        -------
        list
            list of transcripts
        """        
        pass
        # assemble based on LOCUS_TAG?
        # return sorted
    
    def _reset(self):
        """Release memory and reset internal hashes"""
        del self._feature_cache
        gc.collect()
        del gc.garbage[:]
        self._feature_cache = { "mrna_like" : {}, "CDS_like" : {}}
