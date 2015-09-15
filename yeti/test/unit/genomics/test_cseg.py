#!/usr/bin/env python
"""Tests for data structures defined in :py:mod:`yeti.genomics.roitools`
"""
import os
import tempfile
import unittest
import itertools
import copy
from yeti.util.services.mini2to3 import cStringIO

from random import shuffle
from pkg_resources import resource_filename, cleanup_resources
from nose.plugins.attrib import attr

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from yeti.test.ref_files import MINI
from yeti.readers.gff import GTF2_TranscriptAssembler, \
                             GFF3_TranscriptAssembler, \
                             GFF3_Reader
from yeti.readers.bed import BED_Reader

from yeti.genomics.c_roitools import GenomicSegment, positions_to_segments, positionlist_to_segments
from yeti.genomics.c_segmentchain import SegmentChain, Transcript

from yeti.genomics.genome_array import GenomeArray
from yeti.util.io.filters import CommentReader
from yeti.util.services.decorators import skip_if_abstract

import warnings
warnings.simplefilter("ignore",DeprecationWarning)

#===============================================================================
# INDEX: Helper functions & constants
#===============================================================================


ANNOTATION_FILES = {
                    "bed100" : resource_filename("yeti","test/data/annotations/100transcripts.bed"),
                    "gff100" : resource_filename("yeti","test/data/annotations/100transcripts.gff"),
                    "gtf100" : resource_filename("yeti","test/data/annotations/100transcripts.gtf"),
                    "cds100" : resource_filename("yeti","test/data/annotations/100cds.bed"),
                    "utr5_100" : resource_filename("yeti","test/data/annotations/100utr5s.bed"),
                    "utr3_100" : resource_filename("yeti","test/data/annotations/100utr3s.bed"),
                    }

def tearDownModule():
    """Remove test dataset files after unit tests are complete"""
    cleanup_resources()

#===============================================================================
# INDEX: Tests for miscellaneous functions
#===============================================================================

@attr(test="unit")
class TestMiscellaneous(unittest.TestCase):
    """Test suite for miscellaneous methods"""
    def test_positions_to_segments(self):
        # single GenomicSegment, both strands
        self.assertListEqual(positions_to_segments("chrA","+",range(100)),
                             [GenomicSegment("chrA",0,100,"+")])
        self.assertListEqual(positions_to_segments("chrA","-",range(100)),
                             [GenomicSegment("chrA",0,100,"-")])
        self.assertListEqual(positions_to_segments("chrA",".",range(100)),
                             [GenomicSegment("chrA",0,100,".")])
        
        # Multiple intervals, both strands
        self.assertListEqual(positions_to_segments("chrA","+",list(range(100))+list(range(150,200))),
                             [GenomicSegment("chrA",0,100,"+"),GenomicSegment("chrA",150,200,"+")])
        self.assertListEqual(positions_to_segments("chrA","-",list(range(100))+list(range(150,200))),
                             [GenomicSegment("chrA",0,100,"-"),GenomicSegment("chrA",150,200,"-")])
        self.assertListEqual(positions_to_segments("chrA",".",list(range(100))+list(range(150,200))),
                             [GenomicSegment("chrA",0,100,"."),GenomicSegment("chrA",150,200,".")])

        # Multiple intervals, both strands, overlapping positions
        self.assertListEqual(positions_to_segments("chrA","+",list(range(100))+list(range(150,200))+list(range(195,205))),
                             [GenomicSegment("chrA",0,100,"+"),GenomicSegment("chrA",150,205,"+")])
        self.assertListEqual(positions_to_segments("chrA","-",list(range(100))+list(range(150,200))+list(range(195,205))),
                             [GenomicSegment("chrA",0,100,"-"),GenomicSegment("chrA",150,205,"-")])
        self.assertListEqual(positions_to_segments("chrA",".",list(range(100))+list(range(150,200))+list(range(195,205))),
                             [GenomicSegment("chrA",0,100,"."),GenomicSegment("chrA",150,205,".")])
        
        # negative controls
        self.assertNotEqual(positions_to_segments("chrA","+",range(100)),
                            positions_to_segments("chrA","-",range(100)))
        self.assertNotEqual(positions_to_segments("chrA","+",range(100)),
                            positions_to_segments("chrA",".",range(100)))
        self.assertNotEqual(positions_to_segments("chrA",".",range(100)),
                            positions_to_segments("chrA","-",range(100)))

        self.assertNotEqual(positions_to_segments("chrA","+",range(100)),
                            positions_to_segments("chrA","+",range(200)))

        self.assertNotEqual(positions_to_segments("chrA","+",range(100)),
                            positions_to_segments("chrB","+",range(100)))

    def test_positionlist_to_segments(self):
        # single GenomicSegment, both strands
        self.assertListEqual(positionlist_to_segments("chrA","+",list(range(100))),
                             [GenomicSegment("chrA",0,100,"+")])
        self.assertListEqual(positionlist_to_segments("chrA","-",list(range(100))),
                             [GenomicSegment("chrA",0,100,"-")])
        self.assertListEqual(positionlist_to_segments("chrA",".",list(range(100))),
                             [GenomicSegment("chrA",0,100,".")])
        
        # Multiple intervals, both strands
        self.assertListEqual(positionlist_to_segments("chrA","+",list(list(range(100)))+list(range(150,200))),
                             [GenomicSegment("chrA",0,100,"+"),GenomicSegment("chrA",150,200,"+")])
        self.assertListEqual(positionlist_to_segments("chrA","-",list(list(range(100)))+list(range(150,200))),
                             [GenomicSegment("chrA",0,100,"-"),GenomicSegment("chrA",150,200,"-")])
        self.assertListEqual(positionlist_to_segments("chrA",".",list(list(range(100)))+list(range(150,200))),
                             [GenomicSegment("chrA",0,100,"."),GenomicSegment("chrA",150,200,".")])

       
        # negative controls
        self.assertNotEqual(positionlist_to_segments("chrA","+",list(range(100))),
                            positionlist_to_segments("chrA","-",list(range(100))))
        self.assertNotEqual(positionlist_to_segments("chrA","+",list(range(100))),
                            positionlist_to_segments("chrA",".",list(range(100))))
        self.assertNotEqual(positionlist_to_segments("chrA",".",list(range(100))),
                            positionlist_to_segments("chrA","-",list(range(100))))

        self.assertNotEqual(positionlist_to_segments("chrA","+",list(range(100))),
                            positionlist_to_segments("chrA","+",list(range(200))))

        self.assertNotEqual(positionlist_to_segments("chrA","+",list(range(100))),
                            positionlist_to_segments("chrB","+",list(range(100))))

#===============================================================================
# INDEX: Tests for SegmentChain and Transcript
#===============================================================================


@attr(test="unit")
class AbstractSegmentChainHelper(unittest.TestCase):
    """Test suite for :py:class:`yeti.genomics.roitools.SegmentChain`"""

    @staticmethod
    def _to_segmentchain(transcript_list):
        """Convert a list of Transcript objects to SegmentChains"""
        ltmp = []
        for tx in transcript_list:
            ltmp.append(SegmentChain(*tx._segments,ID=tx.get_name()))
        
        assert len(ltmp) == len(transcript_list)
        
        return ltmp
    
    @staticmethod
    def is_identical(ivc1,ivc2):
        """Test for identity between positions of two SegmentChains,
        defining equality as identity of position sets, strand, and
        chromosomes. Attributes are ignored.
        """
        position_test = ivc1.get_position_set() == ivc2.get_position_set()
        strand_test   = ivc1.spanning_segment.strand == ivc2.spanning_segment.strand
        chrom_test    = ivc1.spanning_segment.chrom == ivc2.spanning_segment.chrom
        return len(ivc1) > 0 and len(ivc2) > 0 and position_test and strand_test and chrom_test
    
    @skip_if_abstract
    def test_import_bed_gtf_gff(self):
        """Assert SegmentChains imported from BED, GTF2, GFF3 files are identical"""
        
        # make sure dictionaries have same keysets
        self.assertEquals(set(self.bed_dict.keys()),set(self.gff_dict.keys()),
                          "BED and GFF %s have different names:\nBED    %s\nGFF    %s" % (self.test_class.__name__,
                                                                                          ",".join(self.bed_dict.keys()),
                                                                                          ",".join(self.gff_dict.keys())))
        self.assertEquals(set(self.bed_dict.keys()),set(self.gtf_dict.keys()),
                          "BED and GTF %s have different names:\nBED    %s\nGTF    %s" % (self.test_class.__name__,
                                                                                          ",".join(self.bed_dict.keys()),
                                                                                          ",".join(self.gtf_dict.keys())))
        
        # make sure transcripts with identical names evaluate identically
        for k in set(self.bed_dict.keys()):
            self.assertTrue(self.is_identical(self.bed_dict[k],self.gtf_dict[k]),
                            "%s is different in bed_dict vs gtf_dict: %s,%s" % (k,self.bed_dict[k],self.gtf_dict[k]))
            self.assertTrue(self.is_identical(self.bed_dict[k],self.gff_dict[k]),
                            "%s is different in bed_dict vs gff_dict: %s,%s" % (k,self.bed_dict[k],self.gtf_dict[k]))
        
        # make sure transcripts with non-identical names evaluate non-identically
        c = 0
        for k1, k2 in zip(self.sorted_keys,self.shuffled_keys):
            if k1 != k2:
                c += 1
                err_msg = "Different %s %s and %s evaluate identically." % (self.test_class.__name__,k1,k2)
                self.assertFalse(self.is_identical(self.bed_dict[k1],self.bed_dict[k2]),err_msg)
                self.assertFalse(self.is_identical(self.gtf_dict[k1],self.gtf_dict[k2]),err_msg)
                self.assertFalse(self.is_identical(self.gff_dict[k1],self.gff_dict[k2]),err_msg)
            
        self.assertGreater(c,0)

    @skip_if_abstract    
    def test_export_bed_gtf(self):
        """Test export to BED12 and GTF formats"""
        bed_text = tempfile.NamedTemporaryFile(mode="w",delete=False)
        gtf_text = tempfile.NamedTemporaryFile(mode="w",delete=False)

        for tx in self.bed_dict.values():
            gtf_text.write(tx.as_gtf())
            bed_text.write(tx.as_bed())

        bed_text.close()
        gtf_text.close()

        print("test_export_bed_gtf: able to write")
        bed_transcripts = {X.get_name() : X  for X in BED_Reader(open(bed_text.name),return_type=SegmentChain) } #Transcript) }
        gtf_transcripts = {X.get_name() : X  for X in GTF2_TranscriptAssembler(open(gtf_text.name),
                                                                               add_three_for_stop=False,
                                                                               return_type=SegmentChain) } #Transcript) }

        print("test_export_bed_gtf: able to read back in")
        # Test equality of imported & exported keysets
        self.assertEquals(set(self.bed_dict.keys()),set(bed_transcripts.keys()))
        self.assertEquals(set(self.bed_dict.keys()),set(gtf_transcripts.keys()))

        # Test equality of re-imported exported transcripts
        for txid, tx in self.bed_dict.items():
            err_msg = "Reimported SegmentChain %s does not evaluate identically to reference." % txid
            self.assertTrue(self.is_identical(tx,bed_transcripts[txid]),err_msg)
            self.assertTrue(self.is_identical(tx,gtf_transcripts[txid]),err_msg)

        print("reimported transcripts evaluate equally")
        # Test inequality of those that should be different
        c = 0
        for k1, k2 in zip(self.sorted_keys,self.shuffled_keys):
            if k1 != k2:
                c += 1
                err_msg = "Different SegmentChains %s and %s evaluate identically." % (k1,k2)
                self.assertFalse(self.is_identical(self.bed_dict[k1],bed_transcripts[k2]),err_msg)
                self.assertFalse(self.is_identical(self.bed_dict[k1],gtf_transcripts[k2]),err_msg)
            
        self.assertGreater(c,0)

        print("Testing inequality of items that should be different")
        os.remove(bed_text.name)
        os.remove(gtf_text.name)        
    
    @skip_if_abstract    
    def test_as_gff3_is_empty_for_zerolen(self):
        self.assertEqual("",SegmentChain().as_gff3())
    
    @skip_if_abstract    
    def test_export_all_GTF2_attributes(self):
        inp = """chrI    .    exon    224    790    .    +    .    gene_id "YAL069W"; transcript_id "YAL069W_mRNA"; name "YAL069W"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL069W"; ID "YAL069W_mRNA";
chrI    .    CDS    335    646    .    +    .    gene_id "YAL069W"; transcript_id "YAL069W_mRNA"; name "YAL069W"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL069W"; ID "YAL069W_mRNA";
chrI    .    start_codon    335    338    .    +    .    gene_id "YAL069W"; transcript_id "YAL069W_mRNA"; name "YAL069W"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL069W"; ID "YAL069W_mRNA";
chrI    .    stop_codon    646    649    .    +    .    gene_id "YAL069W"; transcript_id "YAL069W_mRNA"; name "YAL069W"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL069W"; ID "YAL069W_mRNA";
chrI    .    exon    427    933    .    +    .    gene_id "YAL068W-A"; transcript_id "YAL068W-A_mRNA"; name "YAL068W-A"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL068W-A"; ID "YAL068W-A_mRNA";
chrI    .    CDS    538    789    .    +    .    gene_id "YAL068W-A"; transcript_id "YAL068W-A_mRNA"; name "YAL068W-A"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL068W-A"; ID "YAL068W-A_mRNA";
chrI    .    start_codon    538    541    .    +    .    gene_id "YAL068W-A"; transcript_id "YAL068W-A_mRNA"; name "YAL068W-A"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL068W-A"; ID "YAL068W-A_mRNA";
chrI    .    stop_codon    789    792    .    +    .    gene_id "YAL068W-A"; transcript_id "YAL068W-A_mRNA"; name "YAL068W-A"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL068W-A"; ID "YAL068W-A_mRNA";
chrI    .    exon    1666    2280    .    -    .    gene_id "YAL068C"; transcript_id "YAL068C_mRNA"; name "PAU8"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "PAU8%2Cseripauperin PAU8"; ID "YAL068C_mRNA";
chrI    .    CDS    1810    2169    .    -    .    gene_id "YAL068C"; transcript_id "YAL068C_mRNA"; name "PAU8"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "PAU8%2Cseripauperin PAU8"; ID "YAL068C_mRNA";
chrI    .    start_codon    2166    2169    .    -    .    gene_id "YAL068C"; transcript_id "YAL068C_mRNA"; name "PAU8"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "PAU8%2Cseripauperin PAU8"; ID "YAL068C_mRNA";
chrI    .    stop_codon    1807    1810    .    -    .    gene_id "YAL068C"; transcript_id "YAL068C_mRNA"; name "PAU8"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "PAU8%2Cseripauperin PAU8"; ID "YAL068C_mRNA";
chrI    .    exon    2369    2848    .    +    .    gene_id "YAL067W-A"; transcript_id "YAL067W-A_mRNA"; name "YAL067W-A"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL067W-A"; ID "YAL067W-A_mRNA";
chrI    .    CDS    2480    2704    .    +    .    gene_id "YAL067W-A"; transcript_id "YAL067W-A_mRNA"; name "YAL067W-A"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL067W-A"; ID "YAL067W-A_mRNA";
chrI    .    start_codon    2701    2704    .    +    .    gene_id "YAL067W-A"; transcript_id "YAL067W-A_mRNA"; name "YAL067W-A"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL067W-A"; ID "YAL067W-A_mRNA";
chrI    .    stop_codon    2704    2707    .    +    .    gene_id "YAL067W-A"; transcript_id "YAL067W-A_mRNA"; name "YAL067W-A"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "YAL067W-A"; ID "YAL067W-A_mRNA";
chrI    .    exon    7094    9127    .    -    .    gene_id "YAL067C"; transcript_id "YAL067C_mRNA"; name "SEO1"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "SEO1%2Cputative permease SEO1"; ID "YAL067C_mRNA";
chrI    .    CDS    7238    9016    .    -    .    gene_id "YAL067C"; transcript_id "YAL067C_mRNA"; name "SEO1"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "SEO1%2Cputative permease SEO1"; ID "YAL067C_mRNA";
chrI    .    start_codon    9013    9016    .    -    .    gene_id "YAL067C"; transcript_id "YAL067C_mRNA"; name "SEO1"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "SEO1%2Cputative permease SEO1"; ID "YAL067C_mRNA";
chrI    .    stop_codon    7235    7238    .    -    .    gene_id "YAL067C"; transcript_id "YAL067C_mRNA"; name "SEO1"; utr5_source "estimated"; utr3_source "estimated"; gene_aliases "SEO1%2Cputative permease SEO1"; ID "YAL067C_mRNA";""".replace("    ","\t")
        transcripts1 = list(GTF2_TranscriptAssembler(cStringIO.StringIO(inp),return_type=self.test_class,add_three_for_stop=False))
        
        buf = cStringIO.StringIO()
        for tx in transcripts1:
            buf.write(tx.as_gtf())
            # make sure we have all 7 attributes plus the 6 generated upon import
            self.assertEqual(len(tx.attr),13)
        
        buf.seek(0)
        transcripts2 = list(GTF2_TranscriptAssembler(buf,return_type=self.test_class,add_three_for_stop=False))
        self.assertEqual(len(transcripts1),len(transcripts2))
        # on export, we expect to have gene_id, transcript_id, 
        for tx1,tx2 in zip(transcripts1,transcripts2):
            self.assertDictEqual(tx1.attr,tx2.attr)
    
    @skip_if_abstract    
    def test_calculate_GTF2_frame_plus(self):
        # test recalculation of phase on plus strand
        iv1p = GenomicSegment("chrA",100,130,"+")
        iv2p = GenomicSegment("chrA",160,192,"+")
        iv3p = GenomicSegment("chrA",220,252,"+")
        iv4p = GenomicSegment("chrA",280,312,"+")
        iv5p = GenomicSegment("chrA",340,370,"+")
        iv6p = GenomicSegment("chrA",400,401,"+")
        iv7p = GenomicSegment("chrA",430,461,"+")
        iv8p = GenomicSegment("chrA",490,521,"+")
        iv9p = GenomicSegment("chrA",550,580,"+")
        ivc = SegmentChain(iv1p,iv2p,iv3p,iv4p,iv5p,iv6p,iv7p,iv8p,iv9p,type="CDS")
        
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        expected_phases = [0,0,1,2,0,0,2,1,0]
        found_phases = []
        for line in buf.read().strip("\n").split("\n"):
            phase = int(line.split("\t")[7])
            found_phases.append(phase)

        self.assertListEqual(found_phases,expected_phases)

    @skip_if_abstract    
    def test_calculate_GTF2_frame_minus(self):
        # test recalculation of phase on minus strand
        iv1m = GenomicSegment("chrA",100,130,"-")
        iv2m = GenomicSegment("chrA",160,192,"-")
        iv3m = GenomicSegment("chrA",220,252,"-")
        iv4m = GenomicSegment("chrA",280,312,"-")
        iv5m = GenomicSegment("chrA",340,370,"-")
        iv6m = GenomicSegment("chrA",400,401,"-")
        iv7m = GenomicSegment("chrA",430,461,"-")
        iv8m = GenomicSegment("chrA",490,521,"-")
        iv9m = GenomicSegment("chrA",550,580,"-")
        ivc = SegmentChain(iv1m,iv2m,iv3m,iv4m,iv5m,iv6m,iv7m,iv8m,iv9m,type="CDS")
        
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        expected_phases = [0,0,1,2,0,0,2,1,0]
        found_phases = []
        for line in buf.read().strip("\n").split("\n"):
            phase = int(line.split("\t")[7])
            found_phases.append(phase)

        self.assertListEqual(found_phases,expected_phases)

    @skip_if_abstract    
    def test_calculate_GTF2_frame_known_if_multi_exons(self):
        # test recalculation of phase on plus strand
        iv1p = GenomicSegment("chrA",100,130,"+")
        iv2p = GenomicSegment("chrA",160,192,"+")
        iv3p = GenomicSegment("chrA",220,252,"+")
        iv4p = GenomicSegment("chrA",280,312,"+")
        iv5p = GenomicSegment("chrA",340,370,"+")
        iv6p = GenomicSegment("chrA",400,401,"+")
        iv7p = GenomicSegment("chrA",430,461,"+")
        iv8p = GenomicSegment("chrA",490,521,"+")
        iv9p = GenomicSegment("chrA",550,580,"+")
        ivc = SegmentChain(iv1p,iv2p,iv3p,iv4p,iv5p,iv6p,iv7p,iv8p,iv9p,type="CDS",frame="2")
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        expected_phases = [0,0,1,2,0,0,2,1,0]
        found_phases = []
        for line in buf.read().strip("\n").split("\n"):
            phase = int(line.split("\t")[7])
            found_phases.append(phase)

        self.assertListEqual(found_phases,expected_phases)

    @skip_if_abstract    
    def test_export_GTF2_frame_known(self):
        # export 'phase' or 'frame' keys from attr dict if known
        # instead of recalculating
        iv1p = GenomicSegment("chrA",100,130,"+")
        ivc = SegmentChain(iv1p,phase=0,type="CDS")
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        print("test1")
        for line in buf.read().strip("\n").split("\n"):
            self.assertEqual(line.split("\t")[7],"0")

        ivc = SegmentChain(iv1p,phase=1,type="CDS")
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        print("test1")
        for line in buf.read().strip("\n").split("\n"):
            self.assertEqual(line.split("\t")[7],"1")

        ivc = SegmentChain(iv1p,phase=1,type="CDS")
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        print("test1")
        for line in buf.read().strip("\n").split("\n"):
            self.assertEqual(line.split("\t")[7],"1")

        ivc = SegmentChain(iv1p,frame=0,type="CDS")
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        print("test1")
        for line in buf.read().strip("\n").split("\n"):
            self.assertEqual(line.split("\t")[7],"0")

        ivc = SegmentChain(iv1p,frame=1,type="CDS")
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        print("test1")
        for line in buf.read().strip("\n").split("\n"):
            self.assertEqual(line.split("\t")[7],"1")
                
    @skip_if_abstract    
    def test_not_export_GTF2_frame_non_cds(self):
        # make sure "frame" column is set to "." for non-CDS features
        iv1m = GenomicSegment("chrA",100,130,"-")
        iv2m = GenomicSegment("chrA",160,192,"-")
        iv3m = GenomicSegment("chrA",220,252,"-")
        iv4m = GenomicSegment("chrA",280,312,"-")
        iv5m = GenomicSegment("chrA",340,370,"-")
        iv6m = GenomicSegment("chrA",400,401,"-")
        iv7m = GenomicSegment("chrA",430,461,"-")
        iv8m = GenomicSegment("chrA",490,521,"-")
        iv9m = GenomicSegment("chrA",550,580,"-")
        ivc = SegmentChain(iv1m,iv2m,iv3m,iv4m,iv5m,iv6m,iv7m,iv8m,iv9m,type="exon")
        
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        expected_phases = ["."]*9
        found_phases = []
        for line in buf.read().strip("\n").split("\n"):
            phase = line.split("\t")[7]
            found_phases.append(phase)
        
        self.assertListEqual(found_phases,expected_phases)

    @skip_if_abstract    
    def test_not_export_GTF2_frame_known_non_cds(self):
        # make sure "frame" column is set to "." for non-CDS features
        # even if it is defined in keywords
        iv1p = GenomicSegment("chrA",100,130,"+")
        iv2p = GenomicSegment("chrA",160,192,"+")
        ivc = SegmentChain(iv1p,iv2p,phase=0,type="exon")
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        for line in buf.read().strip("\n").split("\n"):
            self.assertEqual(line.split("\t")[7],".")

        ivc = SegmentChain(iv1p,iv2p,frame=0,type="exon")
        buf = cStringIO.StringIO()
        buf.write(ivc.as_gtf())
        buf.seek(0)
        for line in buf.read().strip("\n").split("\n"):
            self.assertEqual(line.split("\t")[7],".")
    
    @skip_if_abstract    
    def test_sort(self):
        """Test sorting of intervals within a SegmentChain
        
        This function is typically automatic and not exposed to users.
        We test it manually here by shuffling intervals within an 
        SegmentChain, asserting that their order be permuted,
        then re-sorting them, and asserting that their original order
        be restored.
        """
        c = 0
        n = 0
        for ivc in self.bed_list:
            if len(ivc) > 1: # can only shuffle if we have >1 interval
                n += 1
                shuffled_ivc = copy.deepcopy(ivc)
                shuffle(shuffled_ivc._segments)
                
                # if we shuffle SegmentChains with 2 intervals,
                # they will 50% of the time retain the same order by chance
                # so we count how many are changed (`c`) and make sure
                # that the total number of passing tests (`n`) exceeds
                # this number
                if shuffled_ivc._segments != ivc._segments:
                    c += 1
                    
                shuffled_ivc.sort()
                self.assertEquals(shuffled_ivc._segments,ivc._segments)
        
        self.assertGreaterEqual(n,c)
    
    @skip_if_abstract    
    def test_len(self):
        """Test reporting of number of intervals"""
        for ivc in self.bed_list:
            self.assertEquals(len(ivc),len(ivc._segments))
    
    @skip_if_abstract    
    def test_identity(self):
        """Test identity of %ss""" % (self.test_class.__name__)
        combos = itertools.product(self.ivcs.keys(),repeat=2)
        for k1, k2 in combos:
            msg = "%s test_identity: '%s' %s and '%s' %s evaluated as %s. Expected %s."
            chain1 = self.ivcs[k1]
            chain2 = self.ivcs[k2]
            print("testing %s (%s), %s (%s)" % (k1,chain1,k2,chain2))
            answer = self.is_identical(chain1,chain2)
            if k1 == k2 and len(chain1) > 0 and len(chain2) > 0:
                self.assertTrue(answer,msg % (self.test_class,k1,chain1,k2,chain2,False,True)  )
            else:
                #self.assertFalse(self.is_identical(self.ivcs[k1],self.ivcs[k2]))
                self.assertFalse(answer,msg % (self.test_class,k1,chain1,k2,chain2,True,False)  )


    def _eq_check(self,test_key,test_func,strand_tests=["sense"]):
        """Test SegmentChains for various stranded equality properties
        
        Parameters
        ----------
        test_key : str
            Key indicating tuples of IVcollections that should pass
            
        test_func : function
            Function that should evaluate to True or False
        
        strand_tests : list<str>
            List that can contain "sense" and/or "antisense"
            to specify the expected relationship between strands
            of the two SegmentChains being compared
        """
        combos = itertools.product(self.ivc_keys,repeat=2)
        chromkeys = ["a","b"]
        strandkeys = ["p","m"]
        
        strand_eq = []
        if "sense" in strand_tests:
            strand_eq += [("p","p"),("m","m")]
        if "antisense" in strand_tests:
            strand_eq += [("p","m"),("m","p")]
            
        strand_eq  = set(strand_eq)
        strand_neq = set(itertools.product(strandkeys,repeat=2)) - strand_eq 
        
        for k1, k2 in combos:
            # empty test should be false
            if k1 == "C" or k2 == "C":
                func = self.assertFalse
                expect = False
            # identities should be true
            elif k1 == k2:
                func = self.assertTrue
                expect = True
            # anything listed in self.tests should be true
            elif (k1,k2) in self.tests[test_key]:
                func = self.assertTrue
                expect = True
            # all else should be false
            else:
                func = self.assertFalse
                expect = False
            for my_chrom in chromkeys:
                # make sure all strand-equal combos evaluate correctly
                for strand1,strand2 in strand_eq:
                    ivc1 = "%s%s%s" % (my_chrom,k1,strand1)
                    ivc2 = "%s%s%s" % (my_chrom,k2,strand2)
                    err_msg = "%s and %s fail `%s` test. Expected: %s" % (ivc1,ivc2,test_key,expect)
                    func(test_func(self.ivcs[ivc1],self.ivcs[ivc2]),err_msg)
                
                # make sure opposite strands are false
                for strand1, strand2 in strand_neq:
                    ivc1 =  "%s%s%s" % (my_chrom,k1,strand1)
                    ivc2 =  "%s%s%s" % (my_chrom,k1,strand2)
                    err_msg = "%s and %s fail opposite-strand `%s` test. Expected: %s" % (ivc1,ivc2,test_key,expect)
                    self.assertFalse(test_func(self.ivcs[ivc1],self.ivcs[ivc2]),err_msg)

            # make sure opposite chromosomes are false
            for strand1,strand2 in strand_eq | strand_neq:
                ivc1 =  "%s%s%s" % ("a",k1,strand1)
                ivc2 =  "%s%s%s" % ("b",k1,strand2)
                err_msg = "%s and %s fail opposite-chromosome `%s` test. Expected: %s" % (ivc1,ivc2,test_key,expect)
                self.assertFalse(test_func(self.ivcs[ivc1],self.ivcs[ivc2]),err_msg)

                ivc1 =  "%s%s%s" % ("b",k1,strand1)
                ivc2 =  "%s%s%s" % ("a",k1,strand2)
                err_msg = "%s and %s fail opposite-chromosome `%s` test. Expected: %s" % (ivc1,ivc2,test_key,expect)
                self.assertFalse(test_func(self.ivcs[ivc1],self.ivcs[ivc2]),err_msg)
                
    @skip_if_abstract    
    def test_covers(self):
        self._eq_check("covers",self.test_class.covers)
    
    @skip_if_abstract    
    def test_shares_segments_with(self):
        self._eq_check("shares_segments_with",self.test_class.shares_segments_with)
    
    @skip_if_abstract    
    def test_overlaps(self):
        self._eq_check("overlaps",self.test_class.overlaps,strand_tests=["sense"])
        
    @skip_if_abstract    
    def test_antisense_overlaps(self):
        self._eq_check("overlaps",self.test_class.antisense_overlaps,strand_tests=["antisense"])

    @skip_if_abstract    
    def test_unstranded_overlaps(self):
        self._eq_check("overlaps",self.test_class.unstranded_overlaps,strand_tests=["sense","antisense"])
    
    @skip_if_abstract    
    def test_contains(self):
        self._eq_check("contains",self.test_class.__contains__)
    
    @skip_if_abstract    
    def test_equals(self):
        self._eq_check("equals",self.test_class.__eq__)
        
    @skip_if_abstract    
    def test_get_antisense(self):
        for ivc in self.bed_list:
            antisense = ivc.get_antisense()
            for iv, asiv in zip(ivc,antisense):
                self.assertEquals(iv.chrom,asiv.chrom)
                self.assertEquals(iv.start,asiv.start)
                self.assertEquals(iv.end,asiv.end)
                self.assertNotEquals(iv.strand,asiv.strand)
    
    @skip_if_abstract    
    def test_get_length(self):
        for ivc in self.bed_list:
            len1 = ivc.get_length()
            len2 = sum(len(X) for X in ivc)
            err_msg = "%s %s fails total length: %s != %s" % (self.test_class.__name__,ivc.get_name(),len1,len2)
            self.assertEquals(len1,len2,err_msg)
    
    @skip_if_abstract    
    def test_add_segments(self):
        ivc = self.test_class()
        
        # add one GenomicSegment
        ivc.add_segments(self.ivs["a1p"])
        self.assertEquals(len(ivc),1)
        self.assertEquals(ivc.get_length(),50)
        
        # add an overlapping GenomicSegment
        ivc.add_segments(self.ivs["a7p"])
        self.assertEquals(len(ivc),1)
        self.assertEquals(ivc.get_length(),75)
        
        # add a non-overlapping GenomicSegment
        ivc.add_segments(self.ivs["a3p"])
        self.assertEquals(len(ivc),2)
        self.assertEquals(ivc.get_length(),175)
        
        # add two GenomicSegments
        ivc.add_segments(self.ivs["a4p"],self.ivs["a5p"])
        self.assertEquals(len(ivc),4)
        self.assertEquals(ivc.get_length(),475)
        
        # add GenomicSegment from the wrong strand
        self.assertRaises(ValueError,
                          ivc.add_segments,self.ivs["a3m"])
        
        # add GenomicSegment from the wrong chromosome
        self.assertRaises(ValueError,
                          ivc.add_segments,self.ivs["b3p"])

        # add GenomicSegment from the wrong chromosome
        self.assertRaises(ValueError,
                          ivc.add_segments,self.ivs["b3m"])
    
    @skip_if_abstract    
    def test_add_masks(self):
        for strand, opp in [("+","-"),("-","+")]:
            ivc = SegmentChain(GenomicSegment("chrA",100,150,strand),
                               GenomicSegment("chrA",250,300,strand))

            mask_a = GenomicSegment("chrA",125,150,strand)
            mask_b = GenomicSegment("chrA",275,300,strand)
            mask_c = GenomicSegment("chrA",275,350,strand)
            
            intron_mask = GenomicSegment("chrA",155,245,strand)

            opp_mask_a = GenomicSegment("chrA",125,150,opp)
            opp_mask_b = GenomicSegment("chrA",275,300,opp)
            opp_mask_c = GenomicSegment("chrA",275,350,opp)

            chrom_mask_a = GenomicSegment("chrB",125,150,strand)
            chrom_mask_b = GenomicSegment("chrB",275,300,strand)
            chrom_mask_c = GenomicSegment("chrB",275,350,strand)
            
            # don't add opposite-strand masks
            self.assertRaises(ValueError,ivc.add_masks,opp_mask_a)
            self.assertEqual(ivc.get_masks(),[])

            # don't add opposite-strand masks
            self.assertRaises(ValueError,ivc.add_masks,opp_mask_a,opp_mask_b,opp_mask_c)
            self.assertEqual(ivc.get_masks(),[])

            # don't add wrong chromosome same-strand masks
            self.assertRaises(ValueError,ivc.add_masks,chrom_mask_a,chrom_mask_b,chrom_mask_c)
            self.assertEqual(ivc.get_masks(),[])
            
            # don't add intron mask
            ivc.add_masks(intron_mask)
            self.assertEqual(ivc.get_masks(),[],"Intron mask added")
            
            # add a single mask
            ivc.add_masks(mask_a)
            self.assertEqual(ivc.get_masks(),[mask_a],"Failure to add single mask")

            # add a single mask again, make sure doesn't change
            ivc.add_masks(mask_a)
            self.assertEqual(ivc.get_masks(),[mask_a],"Failed sequential mask add")
            
            # add a sequential mask, make sure we don't lose the first
            ivc.add_masks(mask_b)
            self.assertEqual(ivc.get_masks(),[mask_a,mask_b],"Failed to add multiple masks")
            
            ivc.reset_masks()
            self.assertEqual(ivc.get_masks(),[])

            # add multiple masks
            ivc.add_masks(mask_a,mask_b)
            self.assertEqual(ivc.get_masks(),[mask_a,mask_b],"Failed to add multiple masks")
            ivc.reset_masks()
            
            # make sure masks are trimmed
            ivc.add_masks(mask_a,mask_c)
            self.assertEqual(ivc.get_masks(),[mask_a,mask_b],"Failed to trim masks")

    @skip_if_abstract    
    def test_get_masks(self):
        for strand in ("+", "-"):
            ivc = SegmentChain(GenomicSegment("chrA",100,150,strand),
                               GenomicSegment("chrA",250,300,strand))
            
            mask_a = GenomicSegment("chrA",125,150,strand)
            mask_b = GenomicSegment("chrA",275,300,strand)
            
            self.assertEqual(ivc.get_masks(),[])
            
            ivc.add_masks(mask_a)
            self.assertEqual(ivc.get_masks(),[mask_a])
            
            ivc.reset_masks()
            self.assertEqual(ivc.get_masks(),[])
            ivc.add_masks(mask_a,mask_b)
            self.assertEqual(ivc.get_masks(),[mask_a,mask_b])
    
    @skip_if_abstract    
    def test_get_masks_as_segmentchain(self):
        for strand in ("+", "-"):
            ivc = SegmentChain(GenomicSegment("chrA",100,150,strand),
                                   GenomicSegment("chrA",250,300,strand))
            
            mask_a = GenomicSegment("chrA",125,150,strand)
            mask_b = GenomicSegment("chrA",275,300,strand)
            
            self.assertEquals(len(ivc.get_masks_as_segmentchain()),0)
            self.assertTrue(isinstance(ivc.get_masks_as_segmentchain(),SegmentChain))
            
            ivc.add_masks(mask_a)
            self.assertEqual(ivc.get_masks_as_segmentchain(),SegmentChain(mask_a))
            
            ivc._mask_segments = []
            ivc.add_masks(mask_a,mask_b)
            self.assertEqual(ivc.get_masks_as_segmentchain(),SegmentChain(mask_a,mask_b))
    
    @skip_if_abstract    
    def test_reset_masks(self):
        for strand in ("+", "-"):
            ivc = SegmentChain(GenomicSegment("chrA",100,150,strand),
                               GenomicSegment("chrA",250,300,strand))

            mask_a = GenomicSegment("chrA",125,150,strand)
            mask_b = GenomicSegment("chrA",275,300,strand)
            
            ivc.add_masks(mask_a,mask_b)
            self.assertEqual(ivc.get_masks(),[mask_a,mask_b])
            ivc.reset_masks()
            self.assertEqual(ivc.get_masks(),[],"Failed to reset masks")

    
    @skip_if_abstract    
    def test_get_junctions(self):
        """Test `get_junctions()`"""
        # Make sure we get nothing from a single-interval segmentchain
        self.assertEquals(self.ivcs["aIp"].get_junctions(),[])
        
        # Make sure we get the right answer for a plus-strand segmentchain
        expected = [GenomicSegment("chrA",200,250,"+")]
        self.assertEquals(expected,self.ivcs["aBp"].get_junctions())

        # Make sure we get the right answer for a minus-strand segmentchain
        expected = [GenomicSegment("chrA",200,250,"-")]
        self.assertEquals(expected,self.ivcs["aBm"].get_junctions())

        # Make sure we get the right answer for a plus-strand segmentchain
        expected = [GenomicSegment("chrA",200,250,"+"),GenomicSegment("chrA",350,500,"+")]
        self.assertEquals(expected,self.ivcs["aAp"].get_junctions())
        
        # Make sure we get the right answer for a minus-strand segmentchain
        expected = [GenomicSegment("chrA",200,250,"-"),GenomicSegment("chrA",350,500,"-")]
        self.assertEquals(expected,self.ivcs["aAm"].get_junctions())
        
    @skip_if_abstract    
    def test_get_segmentchain_coordinate(self):
        """Test conversion of genomic coordinates to %s coordinates""" % (self.test_class.__name__)
        # Test plus strand
        iv1 = GenomicSegment("chrA",20,40,"+")
        iv2 = GenomicSegment("chrA",60,70,"+")
        iv3 = GenomicSegment("chrA",80,90,"+")
        ivc = self.test_class(iv1,iv2,iv3)
        c = 0
        for iv in ivc:
            for x in range(iv.start,iv.end):
                ivc_coordinate = ivc.get_segmentchain_coordinate("chrA",x,"+",stranded=True)
                err_msg = "Plus-stranded segmentchain coordinate incorrect (expected %s, found %s)." % (c,ivc_coordinate)
                self.assertEquals(c,ivc_coordinate,err_msg)
                ivc_coordinate = ivc.get_segmentchain_coordinate("chrA",x,"+",stranded=False)
                err_msg = "Plus-unstranded segmentchain coordinate incorrect (expected %s, found %s)." % (c,ivc_coordinate)
                self.assertEquals(c,ivc_coordinate,err_msg)
                c += 1

        # Test minus strand
        iv1 = GenomicSegment("chrA",20,40,"-")
        iv2 = GenomicSegment("chrA",60,70,"-")
        iv3 = GenomicSegment("chrA",80,90,"-")
        ivc = self.test_class(iv1,iv2,iv3)
        c = 0
        for iv in ivc:
            for x in range(iv.start,iv.end):
                ivc_coordinate = ivc.get_segmentchain_coordinate("chrA",x,"-",stranded=True)
                expected = ivc.get_length()-1-c
                err_msg = "Minus-stranded segmentchain coordinate incorrect (expected %s, found %s)." % (expected,ivc_coordinate)
                self.assertEquals(expected,ivc_coordinate,err_msg)

                ivc_coordinate = ivc.get_segmentchain_coordinate("chrA",x,"-",stranded=False)
                err_msg = "Minus-unstranded segmentchain coordinate incorrect (expected %s, found %s)." % (c,ivc_coordinate)
                self.assertEquals(c,ivc_coordinate,err_msg)
                c += 1    
                    
    @skip_if_abstract    
    def test_get_genomic_coordinate(self):
        #"""Test conversions between genome-centric and SegmentChain-centric coordinates"""
        iv1 = GenomicSegment("chrA",2,3,"+")
        iv2 = GenomicSegment("chrA",15,19,"+") 
        iv3 = GenomicSegment("chrA",20,24,"+")
        iv4 = GenomicSegment("chrA",29,30,"+")
        nv1 = GenomicSegment("chrA",2,3,"-")
        nv2 = GenomicSegment("chrA",15,19,"-")
        nv3 = GenomicSegment("chrA",20,24,"-")
        nv4 = GenomicSegment("chrA",29,30,"-")
        ivca = self.test_class(iv1,iv2,iv3,iv4)
        nvca = self.test_class(nv1,nv2,nv3,nv4)

        # test endpoints
        self.assertEquals(ivca.get_genomic_coordinate(0,stranded=True)[1],iv1.start)
        self.assertEquals(ivca.get_genomic_coordinate(ivca.get_length() - 1,stranded=True)[1],iv4.end - 1)
        self.assertEquals(ivca.get_genomic_coordinate(0,stranded=False)[1],iv1.start)
        self.assertEquals(ivca.get_genomic_coordinate(ivca.get_length() - 1,stranded=False)[1],iv4.end - 1)

        self.assertEquals(nvca.get_genomic_coordinate(0,stranded=True)[1],nv4.end - 1)
        self.assertEquals(nvca.get_genomic_coordinate(nvca.get_length()-1,stranded=True)[1],nv1.start)
        self.assertEquals(nvca.get_genomic_coordinate(0,stranded=False)[1],nv1.start)
        self.assertEquals(nvca.get_genomic_coordinate(nvca.get_length()-1,stranded=False)[1],nv4.end - 1)

        # reconstruct all IVs from individual positions
        for ivc in (ivca,nvca):
            positions = []
            for i in range(ivc.get_length()):
                positions.append(ivc.get_genomic_coordinate(i)[1])
            positions = set(positions)
            new_ivs = positions_to_segments(ivc.chrom,ivc.strand,positions)
            for ivA in new_ivs:
                found = False
                for ivB in ivc:
                    if ivA == ivB:
                        found = True
                        break
                self.assertTrue(found)

        # make sure is inverse function of get_segmentchain_coordinate
        for ivc in (ivca,nvca):
            for i in range(ivc.get_length()):
                x = ivc.get_genomic_coordinate(i,stranded=True)[1]
                self.assertEquals(i,ivc.get_segmentchain_coordinate(ivc.chrom,x,ivc.strand,stranded=True))
                x = ivc.get_genomic_coordinate(i,stranded=False)[1]
                self.assertEquals(i,ivc.get_segmentchain_coordinate(ivc.chrom,x,ivc.strand,stranded=False))
    
    @skip_if_abstract    
    def test_get_subchain(self):
        """Test fetching of subregions of %s as SegmentChains""" % (self.test_class.__name__)
        iv1 = GenomicSegment("chrA",2,3,"+")
        iv2 = GenomicSegment("chrA",15,19,"+") 
        iv3 = GenomicSegment("chrA",20,24,"+")
        iv4 = GenomicSegment("chrA",29,30,"+")
        nv1 = GenomicSegment("chrA",2,3,"-")
        nv2 = GenomicSegment("chrA",15,19,"-")
        nv3 = GenomicSegment("chrA",20,24,"-")
        nv4 = GenomicSegment("chrA",29,30,"-")
        ivca = self.test_class(iv1,iv2,iv3,iv4)
        nvca = self.test_class(nv1,nv2,nv3,nv4)

        # assert that full range reproduces full transcript on plus strand
        whole_range = ivca.get_subchain(0,ivca.get_length())
        self.assertTrue(self.is_identical(whole_range,ivca),"test_get_subchain: expected %s, got %s" % (ivca,whole_range))
        self.assertTrue(ivca.covers(whole_range))

        # assert that full range reproduces full transcript on minus strand
        whole_range = nvca.get_subchain(0,nvca.get_length())
        self.assertTrue(self.is_identical(whole_range,nvca),"test_get_subchain: expected %s, got %s" % (nvca,whole_range))
        self.assertTrue(nvca.covers(whole_range))
        
        # assert they don't cover each other
        self.assertFalse(self.is_identical(ivca,nvca))
        
        # take sub range on plus strand, covering 3 GenomicSegments
        expected_set_plus = set([16, 17, 18, 20, 21, 22, 23, 29])
        sub_plus = ivca.get_subchain(2,27)
        self.assertEquals(sub_plus.get_position_set(),
                          expected_set_plus)
        self.assertEquals(sub_plus.get_position_set(),
                          set(sorted(ivca.get_position_list()[2:27])))
        self.assertTrue(ivca.covers(sub_plus))
        self.assertFalse(sub_plus.covers(ivca))
        
        # take sub range on minus strand, covering 3 GenomicSegments
        expected_set_minus = set([2, 15, 16, 17, 18, 20, 21, 22])
        sub_minus = nvca.get_subchain(2,27)
        self.assertEquals(sub_minus.get_position_set(),
                          expected_set_minus)
        self.assertEquals(sub_minus.get_position_set(),
                          set(sorted(nvca.get_position_list()[::-1][2:27])))
        self.assertTrue(nvca.covers(sub_minus))
        self.assertFalse(sub_plus.covers(nvca))
    
#    @skip_if_abstract    
#    def test_get_counts(self):
#        """Test `get_counts()`, `get_masked_counts()`, and `add_masks()`"""
#        assert False
#        ga = GenomeArray({"chrA":2000})
#        ga[GenomicSegment("chrA",100,200,"+")] = 1
#        ga[GenomicSegment("chrA",250,350,"+")] = 1
#        
#        iv1 = GenomicSegment("chrA",100,150,"+")
#        iv2 = GenomicSegment("chrA",150,200,"+")
#        iv3 = GenomicSegment("chrA",250,350,"+")
#
#        mask = GenomicSegment("chrA",50,125,"+")
#        non_overlap_mask = GenomicSegment("chrA",400,500,"+")
#        
#        ivc1 = self.test_class(iv1,iv2,iv3)
#        ivc2 = self.test_class(iv1,iv2,iv3)
#        
#        pre_unmask_counts = sum(ivc1.get_counts(ga))
#        pre_mask_counts = ivc1.get_masked_counts(ga).sum()
#        self.assertEquals(pre_unmask_counts,200)
#        self.assertEquals(pre_mask_counts,200)
#
#        # add real mask
#        ivc1.add_masks(mask)
#        
#        post_mask_counts = ivc1.get_masked_counts(ga).sum()
#        post_unmask_counts = sum(ivc1.get_counts(ga))
#        self.assertEquals(post_mask_counts,175)
#        self.assertEquals(post_unmask_counts,200)
#        
#        # add non-overlapping mask
#        pre_unmask_counts = sum(ivc2.get_counts(ga))
#        pre_mask_counts = ivc2.get_masked_counts(ga).sum()
#        self.assertEquals(pre_unmask_counts,200)
#        self.assertEquals(pre_mask_counts,200)
#        
#        ivc2.add_masks(non_overlap_mask)
#        post_unmask_counts = sum(ivc2.get_counts(ga))
#        post_mask_counts   = ivc2.get_masked_counts(ga).sum()       
#        self.assertEquals(post_unmask_counts,200)
#        self.assertEquals(post_mask_counts,200)

    @skip_if_abstract    
    def test_masked_total_length(self):
        """Test `get_masked_length()` and `add_masks()`"""
        iv1 = GenomicSegment("chrA",100,150,"+")
        iv2 = GenomicSegment("chrA",150,200,"+")
        iv3 = GenomicSegment("chrA",250,350,"+")

        mask = GenomicSegment("chrA",50,125,"+")
        non_overlap_mask = GenomicSegment("chrA",400,500,"+")
        
        ivc1 = SegmentChain(iv1,iv2,iv3)
        ivc2 = SegmentChain(iv1,iv2,iv3)
        
        pre_unmask_length = ivc1.get_length()
        pre_mask_length = ivc1.get_masked_length()
        self.assertEquals(pre_unmask_length,200)
        self.assertEquals(pre_mask_length,200)

        # add real mask
        ivc1.add_masks(mask)
        
        post_unmask_length = ivc1.get_length()
        post_mask_length = ivc1.get_masked_length()
        self.assertEquals(post_mask_length,175)
        self.assertEquals(post_unmask_length,200)
        
        # add non-overlapping mask
        pre_unmask_length = ivc2.get_length()
        pre_mask_length = ivc2.get_masked_length()
        self.assertEquals(pre_unmask_length,200)
        self.assertEquals(pre_mask_length,200)
        
        ivc2.add_masks(non_overlap_mask)
        post_unmask_length = ivc2.get_length()
        post_mask_length   = ivc2.get_masked_length()
        self.assertEquals(post_unmask_length,200)
        self.assertEquals(post_mask_length,200)

    @skip_if_abstract    
    def test_get_position(self):
        """Test `get_position_set()`, `get_position_list()` and `test_get_masked_position_set()`"""
        iv1 = GenomicSegment("chrA",100,150,"+")
        iv2 = GenomicSegment("chrA",150,200,"+")
        iv3 = GenomicSegment("chrA",250,350,"+")

        mask = GenomicSegment("chrA",50,125,"+")
        non_overlap_mask = GenomicSegment("chrA",400,500,"+")
        
        ivc1 = SegmentChain(iv1,iv2,iv3)
        ivc2 = SegmentChain(iv1,iv2,iv3)
        
        whole_position_list = list(range(100,150)) + list(range(150,200)) + list(range(250,350))
        whole_position_set = set(whole_position_list)
        masked_set = whole_position_set - set(range(50,125))

        pre_unmask_list = ivc1.get_position_list()
        pre_unmask_set  = ivc1.get_position_set()
        pre_unmask_valid_set = ivc1.get_masked_position_set()

        self.assertEquals(pre_unmask_list,whole_position_list)
        self.assertEquals(pre_unmask_set,whole_position_set)
        self.assertEquals(pre_unmask_valid_set,whole_position_set)

        # add real mask
        ivc1.add_masks(mask)
        
        post_unmask_list = ivc1.get_position_list()
        post_unmask_set  = ivc1.get_position_set()
        post_unmask_valid_set = ivc1.get_masked_position_set()

        self.assertEquals(post_unmask_list,whole_position_list)
        self.assertEquals(post_unmask_set,whole_position_set)
        self.assertEquals(post_unmask_valid_set,masked_set)
        
        # add non-overlapping mask
        pre_unmask_list = ivc2.get_position_list()
        pre_unmask_set  = ivc2.get_position_set()
        pre_unmask_valid_set = ivc2.get_masked_position_set()

        self.assertEquals(pre_unmask_list,whole_position_list)
        self.assertEquals(pre_unmask_set,whole_position_set)
        self.assertEquals(pre_unmask_valid_set,whole_position_set)
        
        ivc2.add_masks(non_overlap_mask)
        post_unmask_list = ivc2.get_position_list()
        post_unmask_set  = ivc2.get_position_set()
        post_unmask_valid_set = ivc2.get_masked_position_set()

        self.assertEquals(post_unmask_list,whole_position_list)
        self.assertEquals(post_unmask_set,whole_position_set)
        self.assertEquals(post_unmask_valid_set,whole_position_set)
        
    @skip_if_abstract    
    def test_get_sequence_seqrecord(self):
        """Test `get_sequence()` and `get_fasta()` with a SeqRecord genome"""
        my_seq = "TCTAGA" + 50*"A" + "CCGCGG" + 30*"T"
        genome = { "chrA" : SeqRecord(Seq(my_seq,generic_dna)) }
        
        my_revcomp = str(genome["chrA"].reverse_complement().seq)
        
        iv1p = GenomicSegment("chrA",0,6,"+")
        iv2p = GenomicSegment("chrA",56,62,"+")
        iv3p = GenomicSegment("chrA",0,92,"+")

        iv1m = GenomicSegment("chrA",0,6,"-")
        iv2m = GenomicSegment("chrA",56,62,"-")
        iv3m = GenomicSegment("chrA",0,92,"-")

        ivc1p = self.test_class(iv1p,iv2p,ID="ivc1p")
        ivc1m = self.test_class(iv1m,iv2m,ID="ivc1m")

        self.assertEquals(ivc1p.get_sequence(genome),"TCTAGACCGCGG")
        self.assertEquals(ivc1p.get_fasta(genome),">ivc1p\nTCTAGACCGCGG\n")
        
        self.assertEquals(ivc1m.get_sequence(genome),"CCGCGGTCTAGA")
        self.assertEquals(ivc1m.get_fasta(genome),">ivc1m\nCCGCGGTCTAGA\n")

        ivc2p = self.test_class(iv3p,ID="ivc2p")
        ivc2m = self.test_class(iv3m,ID="ivc2m")
        
        self.assertEquals(ivc2p.get_sequence(genome),my_seq)
        self.assertEquals(ivc2p.get_fasta(genome),">ivc2p\n%s\n" %my_seq)

        self.assertEquals(ivc2m.get_sequence(genome),my_revcomp)
        self.assertEquals(ivc2m.get_fasta(genome),">ivc2m\n%s\n" %my_revcomp)
        
    @skip_if_abstract    
    def test_get_sequence_str(self):
        """Test `get_sequence()` and `get_fasta()` with a string genome"""
        my_seq = "TCTAGA" + 50*"A" + "CCGCGG" + 30*"T"
        genome = { "chrA" : my_seq }
        
        my_revcomp = str(SeqRecord(Seq(genome["chrA"],generic_dna)).reverse_complement().seq)
        
        iv1p = GenomicSegment("chrA",0,6,"+")
        iv2p = GenomicSegment("chrA",56,62,"+")
        iv3p = GenomicSegment("chrA",0,92,"+")

        iv1m = GenomicSegment("chrA",0,6,"-")
        iv2m = GenomicSegment("chrA",56,62,"-")
        iv3m = GenomicSegment("chrA",0,92,"-")

        ivc1p = self.test_class(iv1p,iv2p,ID="ivc1p")
        ivc1m = self.test_class(iv1m,iv2m,ID="ivc1m")

        self.assertEquals(ivc1p.get_sequence(genome),"TCTAGACCGCGG")
        self.assertEquals(ivc1p.get_fasta(genome),">ivc1p\nTCTAGACCGCGG\n")
        
        self.assertEquals(ivc1m.get_sequence(genome),"CCGCGGTCTAGA")
        self.assertEquals(ivc1m.get_fasta(genome),">ivc1m\nCCGCGGTCTAGA\n")

        ivc2p = self.test_class(iv3p,ID="ivc2p")
        ivc2m = self.test_class(iv3m,ID="ivc2m")
        
        self.assertEquals(ivc2p.get_sequence(genome),my_seq)
        self.assertEquals(ivc2p.get_fasta(genome),">ivc2p\n%s\n" %my_seq)

        self.assertEquals(ivc2m.get_sequence(genome),my_revcomp)
        self.assertEquals(ivc2m.get_fasta(genome),">ivc2m\n%s\n" %my_revcomp)

    @skip_if_abstract    
    def test_to_from_str_identity(self):
        """Test import to and from strings"""
        for ivc_id, ivc in self.bed_dict.items():
            err_msg = "%s %s fails string io." % (self.test_class.__name__,ivc_id)
            self.assertTrue(self.is_identical(ivc,self.test_class.from_str(str(ivc))),err_msg)

    @skip_if_abstract    
    def test_to_from_bed_identity(self):
        """Test import to and from BED12 format
        
        NOTE: tests with varying numbers of BED columns (3-12) and with
        or without extra_columns are in yeti.test.unit.readers.test_bed
        """
        for ivc_id, ivc in self.bed_dict.items():
            err_msg = "%s %s fails BED io." % (self.test_class.__name__,ivc_id)
            self.assertTrue(self.is_identical(ivc,self.test_class.from_bed(ivc.as_bed())),err_msg)
 
    @skip_if_abstract    
    def test_from_psl(self):
        ref_transcripts = { X.get_name() : X for X in BED_Reader(open(MINI["bed_file"]),return_type=self.test_class) }
        # PSL file won't have knowledge of CDS, so we set these to none in our reference transcripts
        if self.test_class == Transcript:
            for tx in ref_transcripts.values():
                tx.cds_start = tx.cds_end = None
                tx.cds_genoem_start = tx.cds_genoem_end = None
        for line in open(MINI["psl_file"]):
            if not line.startswith("psLayout") and not line.startswith("match") and not line.startswith("---") and not line.startswith(" "):
                psl = self.test_class.from_psl(line)
                ivc_id = psl.get_name()
                if not "repeat_1" in ivc_id and not "repeat_2" in ivc_id:
                    # repeat regions BLAT to wrong place, for obvious regions. so we ignore them
                    ref = ref_transcripts[ivc_id]
                    err_msg = "SegmentChain %s fails PSL import: \n    bed: %s \n    psl: %s." % (ivc_id, str(ref),str(psl))
                    self.assertTrue(self.is_identical(ref,psl),err_msg)
            
    @skip_if_abstract    
    def test_to_psl_raises_error_if_attributes_not_defined(self):
        for chain in BED_Reader(open(MINI["bed_file"]),return_type=self.test_class):
            self.assertRaises(AttributeError,chain.as_psl)
    
    @skip_if_abstract    
    def test_to_from_psl_identity(self):
        for line in open(MINI["psl_file"]):
            if not line.startswith("psLayout") and not line.startswith("match") and not line.startswith("---") and not line.startswith(" "):
                psl = SegmentChain.from_psl(line)
                line_out = psl.as_psl()
                self.assertEqual(line.strip("\n"),line_out.strip("\n"),"PSL export doesn't match input: \n    %s\n    %s" % (line,line_out)) 
        
    @skip_if_abstract
    def test_as_bed_extra_columns(self):
        attr = { "ID" : "some feature ID",
                     "extra_field_1" : 542,
                     "extra_field_2" : "some extra field",
                   }

        my_chain = self.test_class(GenomicSegment("chrA",100,150,"+"),
                                   GenomicSegment("chrA",500,550,"+"),
                                   **attr)
        self.assertEqual(my_chain.as_bed(extra_columns=["extra_field_1","extra_field_2"]).strip(),
                         "chrA	100	550	some feature ID	0	+	100	100	0,0,0	2	50,50,	0,400,	542	some extra field")

        self.assertEqual(my_chain.as_bed(extra_columns=["extra_field_1","nonexistent_field","extra_field_2"]).strip(),
                         "chrA	100	550	some feature ID	0	+	100	100	0,0,0	2	50,50,	0,400,	542		some extra field")

        my_chain.attr['_bedx_column_order'] = ["extra_field_1","extra_field_2"]
        self.assertEqual(my_chain.as_bed().strip(),
                         "chrA	100	550	some feature ID	0	+	100	100	0,0,0	2	50,50,	0,400,	542	some extra field")
        self.assertEqual(my_chain.as_bed(extra_columns=[]).strip(),
                         "chrA	100	550	some feature ID	0	+	100	100	0,0,0	2	50,50,	0,400,")

#     def test_test_class(self):
#         # make sure all these tests are run in subclasses using e.g. 'Transcript'
#         # instead of inherited 'SegmentChain'
#         print(self.test_class.__name__)
#         assert False

    @skip_if_abstract
    def test_get_position_set(self):
        chain1 = self.test_class(GenomicSegment("chrA",50,100,"+"),
                                 GenomicSegment("chrA",120,130,"+"))
        self.assertEqual(chain1.get_position_set(), set(range(50,100)) | set(range(120,130)),
                         "%s failed to get correct position set on plus strand." % self.test_class.__name__)

        chain2 = self.test_class(GenomicSegment("chrA",50,100,"-"),
                                 GenomicSegment("chrA",120,130,"-"))
        self.assertEqual(chain1.get_position_set(), set(range(50,100)) | set(range(120,130)),
                         "%s failed to get correct position set on minus strand." % self.test_class.__name__)

    @skip_if_abstract
    def test_get_position_list(self):
        chain1 = self.test_class(GenomicSegment("chrA",50,100,"+"),
                                 GenomicSegment("chrA",120,130,"+"))
        self.assertEqual(chain1.get_position_list(), list(range(50,100)) + list(range(120,130)),
                         "%s failed to get correct position list on plus strand." % self.test_class.__name__)

        chain2 = self.test_class(GenomicSegment("chrA",50,100,"-"),
                                 GenomicSegment("chrA",120,130,"-"))
        self.assertEqual(chain1.get_position_list(), list(range(50,100)) + list(range(120,130)),
                         "%s failed to get correct position list on minus strand." % self.test_class.__name__)




@attr(test="unit")
class TestSegmentChain(AbstractSegmentChainHelper):
    """Test suite for :py:class:`yeti.genomics.roitools.SegmentChain`"""

    @classmethod
    def setUpClass(cls):
        cls.test_class = SegmentChain
        
        cls.bed_list = list(BED_Reader(CommentReader(open(ANNOTATION_FILES["bed100"])),return_type=SegmentChain))
        cls.gff_list = list(GFF3_TranscriptAssembler(open(ANNOTATION_FILES["gff100"]),
                                                   exon_types=["exon"],add_three_for_stop=False,
                                                   return_type=SegmentChain))
        cls.gtf_list = list(GTF2_TranscriptAssembler(open(ANNOTATION_FILES["gtf100"]),
                                                   return_type=SegmentChain,
                                                   add_three_for_stop=False))      
        
        cls.gff_list = cls._to_segmentchain(cls.gff_list)
        cls.gtf_list = cls._to_segmentchain(cls.gff_list)
        
        cls.bed_dict = { X.get_name() : X for X in cls.bed_list }
        cls.gff_dict = { X.get_name() : X for X in cls.gff_list }
        cls.gtf_dict = { X.get_name() : X for X in cls.gtf_list }
        
        cls.sorted_keys = sorted(cls.bed_dict.keys())
        cls.shuffled_keys = cls.sorted_keys[:]
        shuffle(cls.shuffled_keys)
        
        cls.ivs  = {}
        cls.ivcs = {}

        cls.ivs["a1p"] = GenomicSegment("chrA",100,150,"+")
        cls.ivs["a2p"] = GenomicSegment("chrA",150,200,"+")
        cls.ivs["a3p"] = GenomicSegment("chrA",250,350,"+")
        cls.ivs["a4p"] = GenomicSegment("chrA",500,700,"+")
        cls.ivs["a5p"] = GenomicSegment("chrA",720,820,"+")
        cls.ivs["a6p"] = GenomicSegment("chrA",100,150,"+")
        cls.ivs["a7p"] = GenomicSegment("chrA",75,125, "+")

        cls.ivs["a1m"] = GenomicSegment("chrA",100,150,"-")
        cls.ivs["a2m"] = GenomicSegment("chrA",150,200,"-")
        cls.ivs["a3m"] = GenomicSegment("chrA",250,350,"-")
        cls.ivs["a4m"] = GenomicSegment("chrA",500,700,"-")
        cls.ivs["a5m"] = GenomicSegment("chrA",720,820,"-")
        cls.ivs["a6m"] = GenomicSegment("chrA",100,150,"-")
        cls.ivs["a7m"] = GenomicSegment("chrA",75,125, "-")
        
        cls.ivs["b1p"] = GenomicSegment("chrB",100,150,"+")
        cls.ivs["b2p"] = GenomicSegment("chrB",150,200,"+")
        cls.ivs["b3p"] = GenomicSegment("chrB",250,350,"+")
        cls.ivs["b4p"] = GenomicSegment("chrB",500,700,"+")
        cls.ivs["b5p"] = GenomicSegment("chrB",720,820,"+")
        cls.ivs["b6p"] = GenomicSegment("chrB",100,150,"+")
        cls.ivs["b7p"] = GenomicSegment("chrB",75,125, "+")

        cls.ivs["b1m"] = GenomicSegment("chrB",100,150,"-")
        cls.ivs["b2m"] = GenomicSegment("chrB",150,200,"-")
        cls.ivs["b3m"] = GenomicSegment("chrB",250,350,"-")
        cls.ivs["b4m"] = GenomicSegment("chrB",500,700,"-")
        cls.ivs["b5m"] = GenomicSegment("chrB",720,820,"-")
        cls.ivs["b6m"] = GenomicSegment("chrB",100,150,"-")
        cls.ivs["b7m"] = GenomicSegment("chrB",75,125, "-")
        
        ivs = cls.ivs
        cls.ivcs["aAp"] = SegmentChain(ivs["a1p"],ivs["a2p"],ivs["a3p"],ivs["a4p"])
        cls.ivcs["aBp"] = SegmentChain(ivs["a2p"],ivs["a3p"])
        cls.ivcs["aCp"] = SegmentChain()
        cls.ivcs["aDp"] = SegmentChain(ivs["a2p"],ivs["a3p"],ivs["a4p"],ivs["a5p"])
        cls.ivcs["aHp"] = SegmentChain(ivs["a5p"])
        cls.ivcs["aIp"] = SegmentChain(ivs["a7p"])

        cls.ivcs["aAm"] = SegmentChain(ivs["a1m"],ivs["a2m"],ivs["a3m"],ivs["a4m"])
        cls.ivcs["aBm"] = SegmentChain(ivs["a2m"],ivs["a3m"])
        cls.ivcs["aCm"] = SegmentChain()
        cls.ivcs["aDm"] = SegmentChain(ivs["a2m"],ivs["a3m"],ivs["a4m"],ivs["a5m"])
        cls.ivcs["aHm"] = SegmentChain(ivs["a5m"])
        cls.ivcs["aIm"] = SegmentChain(ivs["a7m"])

        cls.ivcs["bAp"] = SegmentChain(ivs["b1p"],ivs["b2p"],ivs["b3p"],ivs["b4p"])
        cls.ivcs["bBp"] = SegmentChain(ivs["b2p"],ivs["b3p"])
        cls.ivcs["bCp"] = SegmentChain()
        cls.ivcs["bDp"] = SegmentChain(ivs["b2p"],ivs["b3p"],ivs["b4p"],ivs["b5p"])
        cls.ivcs["bHp"] = SegmentChain(ivs["b5p"])
        cls.ivcs["bIp"] = SegmentChain(ivs["b7p"])

        cls.ivcs["bAm"] = SegmentChain(ivs["b1m"],ivs["b2m"],ivs["b3m"],ivs["b4m"])
        cls.ivcs["bBm"] = SegmentChain(ivs["b2m"],ivs["b3m"])
        cls.ivcs["bCm"] = SegmentChain()
        cls.ivcs["bDm"] = SegmentChain(ivs["b2m"],ivs["b3m"],ivs["b4m"],ivs["b5m"])
        cls.ivcs["bHm"] = SegmentChain(ivs["b5m"])
        cls.ivcs["bIm"] = SegmentChain(ivs["b7m"])

        cls.tests = {}
        cls.tests["overlaps"] = [("A","B"),
                                  ("B","A"),
                                  ("A","D"),
                                  ("D","A"),
                                  ("B","D"),
                                  ("D","B"),
                                  ("H","D"),
                                  ("D","H"),
                                  ("I","A"),
                                  ("A","I")
                                 ]

        cls.tests["covers"] = [("A","B"),
                                ("D","B"),
                                ("D","H"),
                               ]

        cls.tests["shares_segments_with"] = [("A","B"),
                                         ("B","A"),
                                         ("A","D"),
                                         ("D","A"),
                                         ("B","D"),
                                         ("D","B"),
                                         ("D","H"),
                                         ("H","D"),
                                        ]
        
        cls.tests["contains"] = [("A","B"),
                                 ("D","B"),
                                 ("D","H"),
                                # ("D","A")
                                 ]
        
        cls.tests["equals"] = [("A","A"),
                                ("B","B"),
                                ("D","D"),
                                ("H","H"),
                                ("I","I"),
                                ]
        
        cls.ivc_keys = list("ABCDHI")

    def test_as_gff3(self):
        inp = """chrI    SGD    intron    87388    87500    .    +    .    Parent=YAL030W_mRNA;Name=YAL030W_intron;orf_classification=Verified
chrII    SGD    five_prime_UTR_intron    45645    45977    .    +    .    Parent=YBL092W_mRNA;Name=YBL092W_five_prime_UTR_intron;orf_classification=Verified
chrII    SGD    intron    142750    142846    .    -    .    Parent=YBL040C_mRNA;Name=YBL040C_intron;orf_classification=Verified
chrII    SGD    intron    186353    186427    .    -    .    Parent=YBL018C_mRNA;Name=YBL018C_intron;orf_classification=Verified
chrII    SGD    transposable_element_gene    259869    265140    .    +    .    ID=YBR012W-B;Name=YBR012W-B;Alias=gag-pol fusion protein;Ontology_term=GO:0000943,GO:0003723,GO:0003887,GO:0003964,GO:0004540,GO:0005634,GO:0008233,GO:0032197;Note=Retrotransposon TYA Gag and TYB Pol genes%3B transcribed/translated as one unit%3B polyprotein is processed to make a nucleocapsid-like protein (Gag)%2C reverse transcriptase (RT)%2C protease (PR)%2C and integrase (IN)%3B similar to retroviral genes;display=Retrotransposon TYA Gag and TYB Pol genes;dbxref=SGD:S000002155
chrII    SGD    plus_1_translational_frameshift    261174    261174    .    +    .    Parent=YBR012W-B;Name=YBR012W-B_plus_1_translational_frameshift;dbxref=SGD:S000002155
chrII    SGD    intron    407028    407122    .    -    .    Parent=YBR082C_mRNA;Name=YBR082C_intron;orf_classification=Verified
chrII    SGD    intron    606281    606668    .    +    .    Parent=YBR191W_mRNA;Name=YBR191W_intron;orf_classification=Verified
chrII    SGD    intron    726918    727011    .    -    .    Parent=YBR255C-A_mRNA;Name=YBR255C-A_intron;orf_classification=Uncharacterized
chrIII    SGD    intron    177907    178213    .    -    .    Parent=YCR031C_mRNA;Name=YCR031C_intron;orf_classification=Verified
chrIV    SGD    intron    254975    255044    .    -    .    Parent=YDL115C_mRNA;Name=YDL115C_intron;orf_classification=Verified
""".replace("    ","\t")

        buf1 = cStringIO.StringIO()
        buf2 = cStringIO.StringIO()
        features1 = list(GFF3_Reader(cStringIO.StringIO(inp)))
        for feature in features1:
            buf1.write(feature.as_gff3())
            buf2.write(feature.as_gff3())
        
        # test import/export equivalence on feature side
        buf1.seek(0)
        features2 = []
        for feature in GFF3_Reader(buf1):
            features2.append(feature)
        for f1, f2 in zip(features1,features2):
            self.assertTrue(self.is_identical(f1,f2))
            self.assertSetEqual(set(f1.attr.keys()),set(f2.attr.keys()))
            for k in f1.attr.keys():
                self.assertEqual(f1.attr[k],f2.attr[k])
        
        # test import/export equivalaence on text side
        buf2.seek(0)
        inp2 = buf2.getvalue()
        lines1 = inp.split("\n")
        lines2 = inp2.split("\n")
        for l1, l2 in zip(lines1,lines2):
            items1 = l1.split("\t")
            items2 = l2.split("\t")
            self.assertListEqual(items1[:8],items2[:8])
            self.assertSetEqual(set(items1[-1].strip(";").split(";")),set(items2[-1].strip(";").split(";")))
    
    def test_as_gff3_on_multiivc_raises_exception(self):
        ivc = self.ivcs["aAp"]
        self.assertRaises(AttributeError,ivc.as_gff3)
        

