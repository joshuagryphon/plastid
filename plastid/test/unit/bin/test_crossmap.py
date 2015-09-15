#!/usr/bin/env python
"""Tests for methods in :py:mod:`plastid.bin.crossmap`
"""
import unittest
from plastid.util.services.mini2to3 import cStringIO
from nose.plugins.attrib import attr
from Bio import SeqIO
from plastid.genomics.roitools import GenomicSegment, SegmentChain
from plastid.util.services.exceptions import MalformedFileError
from plastid.bin.crossmap import simulate_reads, \
                                     FastaNameReader, \
                                     revcomp_mask_chain, \
                                     fa_to_bed

#===============================================================================
# INDEX: test case
#===============================================================================

@attr(test="unit")
class TestCrossmap(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.ivcs = {
                     "plus" : [SegmentChain(GenomicSegment("chrA",0,100,"+")),
                               SegmentChain(GenomicSegment("chrA",50,100,"+")),
                               SegmentChain(GenomicSegment("chrA",50,51,"+"))
                              ],
                      "minus_k25_off0" : [SegmentChain(GenomicSegment("chrA",0+25-1 ,100+25-1,"-")),
                                          SegmentChain(GenomicSegment("chrA",50+25-1,100+25-1,"-")),
                                          SegmentChain(GenomicSegment("chrA",50+25-1,51+25-1 ,"-"))
                                         ],
                      "minus_k50_off0" : [SegmentChain(GenomicSegment("chrA",0+50-1 ,100+50-1,"-")),
                                          SegmentChain(GenomicSegment("chrA",50+50-1,100+50-1,"-")),
                                          SegmentChain(GenomicSegment("chrA",50+50-1,51+50-1 ,"-"))
                                         ],
                      "minus_k25_off10" : [SegmentChain(GenomicSegment("chrA",0+25-1-2*10 ,100+25-1-2*10,"-")),
                                           SegmentChain(GenomicSegment("chrA",50+25-1-2*10,100+25-1-2*10,"-")),
                                           SegmentChain(GenomicSegment("chrA",50+25-1-2*10,51+25-1-2*10 ,"-"))
                                          ],
                      "minus_k50_off10" : [SegmentChain(GenomicSegment("chrA",0+50-1-2*10 ,100+50-1-2*10,"-")),
                                           SegmentChain(GenomicSegment("chrA",50+50-1-2*10,100+50-1-2*10,"-")),
                                           SegmentChain(GenomicSegment("chrA",50+50-1-2*10,51+50-1-2*10 ,"-"))
                                          ],
                    }

    def test_simulate_reads(self):  
        fh = cStringIO.StringIO()
        genome = SeqIO.parse(cStringIO.StringIO(SHORT_FASTA),"fasta")
        for seq in genome:
            simulate_reads(seq,k=25,fh=fh)
        fh.seek(0)
        reads = fh.read()
        self.assertEqual(reads,SHORT_FASTA_KMERS)
    
    def test_fasta_name_reader(self):
        # make sure we get out only names of sequences from FASTA file
        reader = FastaNameReader(cStringIO.StringIO(SHORT_FASTA))
        self.assertEqual(next(reader),"chr50a")
        self.assertEqual(next(reader),"chr30b")
        self.assertRaises(StopIteration,reader.next)
    
    def test_revcomp_mask_chain_no_offset(self):
        for k in (25,50):
            for offset in (0,10):
                revlist = [revcomp_mask_chain(X,k,offset) for X in self.ivcs["plus"]]
                self.assertListEqual(self.ivcs["minus_k%s_off%s" % (k,offset)],revlist)
    
    def test_fa_to_bed(self):
        # test with and without offsets
        # start with block of FASTA formatted sequence:
        #     two blocks on same chromosome
        #     additional block on next chromosome
        #     additional block on next chromosome
        # two blocks on chrA, one on chrB
        # test without offset
        self.maxDiff = None

        reader = cStringIO.StringIO(CROSSMAP_BLOCK)

        # test without offset
        blocks1 = list(fa_to_bed(reader,25,offset=0))
        self.assertListEqual(blocks1,CROSSMAP1)

        # test with offset
        reader = cStringIO.StringIO(CROSSMAP_BLOCK)
        blocks2 = list(fa_to_bed(reader,25,offset=1000))
        self.assertListEqual(blocks2,CROSSMAP2)
        
    def test_fa_to_bed_throws_expected_error(self):
        # test with and without offsets
        # start with unsorted FASTA block
        
        # grab chrA only and randomize order
        reads = cStringIO.StringIO(SHORT_FASTA_KMERS).read().split("\n")
        reads = reads[40:50] + reads[30:40] + reads[0:16] + reads[20:30]
        reads = "\n".join(reads)
        reader = cStringIO.StringIO(reads)

        # fa_to_bed is a generator, so we need to create a callable to make
        # it a list for it to actually raise the error
        tfunc = lambda x,y,z: list(fa_to_bed(x,y,offset=z))
         
        self.assertRaises(MalformedFileError,tfunc,reader,25,0)

        reader = cStringIO.StringIO(reads)
        self.assertRaises(MalformedFileError,tfunc,reader,25,1000)
        reader  = cStringIO.StringIO(reads)
        self.assertRaises(MalformedFileError,tfunc,reader,15,0)


#===============================================================================
# INDEX: test data
#===============================================================================

SHORT_FASTA = """>chr50a
CCGTGATATGACCTAGGTCGAGAGCTAAGCCTCAATGATGCGCTGGCGAT
>chr30b
CCCTCCTTCCGCTGGCCCCGACTGCCCCAG"""

SHORT_FASTA_KMERS=""">chr50a:0(+)
CCGTGATATGACCTAGGTCGAGAGC
>chr50a:1(+)
CGTGATATGACCTAGGTCGAGAGCT
>chr50a:2(+)
GTGATATGACCTAGGTCGAGAGCTA
>chr50a:3(+)
TGATATGACCTAGGTCGAGAGCTAA
>chr50a:4(+)
GATATGACCTAGGTCGAGAGCTAAG
>chr50a:5(+)
ATATGACCTAGGTCGAGAGCTAAGC
>chr50a:6(+)
TATGACCTAGGTCGAGAGCTAAGCC
>chr50a:7(+)
ATGACCTAGGTCGAGAGCTAAGCCT
>chr50a:8(+)
TGACCTAGGTCGAGAGCTAAGCCTC
>chr50a:9(+)
GACCTAGGTCGAGAGCTAAGCCTCA
>chr50a:10(+)
ACCTAGGTCGAGAGCTAAGCCTCAA
>chr50a:11(+)
CCTAGGTCGAGAGCTAAGCCTCAAT
>chr50a:12(+)
CTAGGTCGAGAGCTAAGCCTCAATG
>chr50a:13(+)
TAGGTCGAGAGCTAAGCCTCAATGA
>chr50a:14(+)
AGGTCGAGAGCTAAGCCTCAATGAT
>chr50a:15(+)
GGTCGAGAGCTAAGCCTCAATGATG
>chr50a:16(+)
GTCGAGAGCTAAGCCTCAATGATGC
>chr50a:17(+)
TCGAGAGCTAAGCCTCAATGATGCG
>chr50a:18(+)
CGAGAGCTAAGCCTCAATGATGCGC
>chr50a:19(+)
GAGAGCTAAGCCTCAATGATGCGCT
>chr50a:20(+)
AGAGCTAAGCCTCAATGATGCGCTG
>chr50a:21(+)
GAGCTAAGCCTCAATGATGCGCTGG
>chr50a:22(+)
AGCTAAGCCTCAATGATGCGCTGGC
>chr50a:23(+)
GCTAAGCCTCAATGATGCGCTGGCG
>chr50a:24(+)
CTAAGCCTCAATGATGCGCTGGCGA
>chr50a:25(+)
TAAGCCTCAATGATGCGCTGGCGAT
>chr30b:0(+)
CCCTCCTTCCGCTGGCCCCGACTGC
>chr30b:1(+)
CCTCCTTCCGCTGGCCCCGACTGCC
>chr30b:2(+)
CTCCTTCCGCTGGCCCCGACTGCCC
>chr30b:3(+)
TCCTTCCGCTGGCCCCGACTGCCCC
>chr30b:4(+)
CCTTCCGCTGGCCCCGACTGCCCCA
>chr30b:5(+)
CTTCCGCTGGCCCCGACTGCCCCAG
"""

CROSSMAP1 = [
    (
        SegmentChain(GenomicSegment("chr50a",1,10,"+")),
        SegmentChain(GenomicSegment("chr50a",1+25-1,10+25-1,"-")),
    ),
    (
        SegmentChain(GenomicSegment("chr50a",19,26,"+")),
        SegmentChain(GenomicSegment("chr50a",19+25-1,26+25-1,"-")),
    ),
    (
        SegmentChain(GenomicSegment("chr30b",0,6,"+")),
        SegmentChain(GenomicSegment("chr30b",0+25-1,6+25-1,"-")),
    )
]

CROSSMAP2 = [
    (
        SegmentChain(GenomicSegment("chr50a",1+1000,10+1000,"+")),
        SegmentChain(GenomicSegment("chr50a",1+25-1-1000,10+25-1-1000,"-")),
    ),
    (
        SegmentChain(GenomicSegment("chr50a",19+1000,26+1000,"+")),
        SegmentChain(GenomicSegment("chr50a",19+25-1-1000,26+25-1-1000,"-")),
    ),
    (
        SegmentChain(GenomicSegment("chr30b",0+1000,6+1000,"+")),
        SegmentChain(GenomicSegment("chr30b",0+25-1-1000,6+25-1-1000,"-")),
    )             
]

# two blocks on chrA, one beginning at position 1, one beginning internally
# one block on chrB, beginning at position 0
CROSSMAP_BLOCK=""">chr50a:1(+)
CGTGATATGACCTAGGTCGAGAGCT
>chr50a:2(+)
GTGATATGACCTAGGTCGAGAGCTA
>chr50a:3(+)
TGATATGACCTAGGTCGAGAGCTAA
>chr50a:4(+)
GATATGACCTAGGTCGAGAGCTAAG
>chr50a:5(+)
ATATGACCTAGGTCGAGAGCTAAGC
>chr50a:6(+)
TATGACCTAGGTCGAGAGCTAAGCC
>chr50a:7(+)
ATGACCTAGGTCGAGAGCTAAGCCT
>chr50a:8(+)
TGACCTAGGTCGAGAGCTAAGCCTC
>chr50a:9(+)
GACCTAGGTCGAGAGCTAAGCCTCA
>chr50a:19(+)
GAGAGCTAAGCCTCAATGATGCGCT
>chr50a:20(+)
AGAGCTAAGCCTCAATGATGCGCTG
>chr50a:21(+)
GAGCTAAGCCTCAATGATGCGCTGG
>chr50a:22(+)
AGCTAAGCCTCAATGATGCGCTGGC
>chr50a:23(+)
GCTAAGCCTCAATGATGCGCTGGCG
>chr50a:24(+)
CTAAGCCTCAATGATGCGCTGGCGA
>chr50a:25(+)
TAAGCCTCAATGATGCGCTGGCGAT
>chr30b:0(+)
CCCTCCTTCCGCTGGCCCCGACTGC
>chr30b:1(+)
CCTCCTTCCGCTGGCCCCGACTGCC
>chr30b:2(+)
CTCCTTCCGCTGGCCCCGACTGCCC
>chr30b:3(+)
TCCTTCCGCTGGCCCCGACTGCCCC
>chr30b:4(+)
CCTTCCGCTGGCCCCGACTGCCCCA
>chr30b:5(+)
CTTCCGCTGGCCCCGACTGCCCCAG
"""    
