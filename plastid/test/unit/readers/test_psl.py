#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.readers.blat`

References
----------
http://pombe.nci.nih.gov/genome/goldenPath/help/blatSpec.html
"""
import unittest
import numpy
from plastid.util.services.mini2to3 import cStringIO
from nose.plugins.attrib import attr
from plastid.genomics.roitools import SegmentChain
from plastid.readers.psl import PSL_Reader, BundledPSL_Reader
from plastid.readers.bed import BED_Reader
from plastid.test.ref_files import MINI

@attr(test="unit")
class TestPSL_Reader(unittest.TestCase):
     
    def test_iter_finds_all(self):
        bed_transcripts = [X for X in BED_Reader(open(MINI["bed_file"]),return_type=SegmentChain) if "intron" not in X.get_name()]
        bed_transcripts = [X for X in bed_transcripts if "repeat" not in X.get_name()]
        psl_transcripts = [X for X in PSL_Reader(open(MINI["psl_file"])) if "repeat" not in X.get_name()]
        
        self.assertEqual(len(bed_transcripts),len(psl_transcripts),"Length mismatch: %s vs %s" %(len(bed_transcripts),len(psl_transcripts)))
        self.assertGreater(len(bed_transcripts),0)

    def test_iter_reads_correct(self):
        bed_transcripts = [X for X in BED_Reader(open(MINI["bed_file"]),return_type=SegmentChain) if "intron" not in X.get_name()]
        bed_transcripts = [X for X in bed_transcripts if "repeat" not in X.get_name()]
        psl_transcripts = [X for X in PSL_Reader(open(MINI["psl_file"])) if "repeat" not in X.get_name()]        
        for bed, psl in zip(bed_transcripts,psl_transcripts):
            bed_name = bed.get_name()
            psl_name = psl.get_name()
            self.assertEqual(bed_name,psl_name,"Name mismatch: %s vs %s" % (bed_name,psl_name))
            psl_positions = psl.get_position_set()
            bed_positions = bed.get_position_set()
            self.assertEqual(bed_positions,psl_positions,"Position mismatch: %s vs %s" % (bed_positions,psl_positions))
            bed_chrom = bed.spanning_segment.chrom
            psl_chrom = psl.spanning_segment.chrom
            self.assertEqual(bed_chrom,psl_chrom,"Chromosome mismatch: %s vs %s" % (bed_chrom,psl_chrom) )
            bed_strand = bed.spanning_segment.strand
            psl_strand = psl.spanning_segment.strand
            self.assertEqual(bed_strand,psl_strand,"Strand mismatch: %s vs %s" % (bed_strand,psl_strand) )


@attr(test="unit")
class TestBundledPSL_Reader(unittest.TestCase):
    
    def test_iter_bundles(self):
        # create a list of entries with same query name, so that these will be bundled by reader
        psl_transcripts = list(PSL_Reader(open(MINI["psl_file"])))
        repeats = numpy.random.randint(1,high=10,size=len(psl_transcripts))
        stmp = ""
        for tx, r in zip(psl_transcripts,repeats):
            for _ in range(r):
                stmp += tx.as_psl()
            
        fakefile = cStringIO.StringIO(stmp)
        reader = BundledPSL_Reader(fakefile)
        for n,tx_group in enumerate(reader):
            # make sure groups are correct length
            self.assertEqual(len(tx_group),repeats[n])
            
            # make sure each group has the correct query name
            names = [X.get_name() for X in tx_group]
            self.assertEqual(len(set(names)),1)
            self.assertEqual(names[0],psl_transcripts[n].get_name())
        
        self.assertEqual(n+1,len(repeats))
        self.assertEqual(n+1,len(psl_transcripts))
