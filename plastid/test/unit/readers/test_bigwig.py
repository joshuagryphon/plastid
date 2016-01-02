import os
from nose.tools import assert_less_equal, assert_raises, assert_dict_equal
from pkg_resources import resource_filename
from plastid.readers.bigwig import BigWigReader
from plastid.genomics.roitools import GenomicSegment
from plastid.genomics.genome_array import GenomeArray
from plastid.util.services.decorators import skip_if_abstract


TOL=1e-5

base_path = resource_filename("plastid","test/data")
wigfile = os.path.join(base_path,"command_line","gen_reads_center_12_fw.wig")
bigwigfile = os.path.join(base_path,"command_line","gen_reads_center_12_fw.bw")


class AbstractTestBBIFile():

    @classmethod
    def setUpClass(cls):
        cls.bw = BigWigReader(bigwigfile) 

    @skip_if_abstract        
    def test_chrom_sizes(self):
        found = self.bw.chroms
        expected = {}
        for line in open(os.path.join(base_path,"annotations","sacCer3.sizes")):
            k,v = line.strip().split("\t")
            expected[k] = int(v)
        
        # these two happen not to be in the dataset
        expected.pop("chrVI")
        expected.pop("chrmt")
            
        assert_dict_equal(expected,found)

    @skip_if_abstract
    def test_no_crash_if_file_not_exist(self):
        with assert_raises(IOError) as _:
            _ = self.reader_class("non_existant_file")


class TestBigWigReader(AbstractTestBBIFile):
    
    @classmethod
    def setUpClass(cls):
        cls.bw = BigWigReader(bigwigfile)
        cls.reader_class = BigWigReader 
    
    def check_vals_against_wig(self,expected,found):
        diff = abs(expected-found)
        maxdiff = diff.max()
        maxloc = diff.argmax()
        msg = "Maximum difference found between BigWig and Wiggle (%s) is at position %s and exceeded tolerance (%s).\n" % (maxdiff,maxloc,TOL)
        msg += "At that position, expected %s, got %s." % (expected[maxloc],found[maxloc])
        assert_less_equal(maxdiff,TOL,msg)
        
    # NOTE: this test relies on WiggleReader being correct
    def test_vals_against_wig(self):
        ga = GenomeArray()
        chrdict = {
            'chrI': 230218,
            'chrII': 813184,
            'chrIII': 316620,
            'chrIV': 1531933,
            'chrIX': 439888,
            'chrV': 576874,
            'chrVII': 1090940,
            'chrVIII': 562643,
            'chrX': 745751,
            'chrXI': 666816,
            'chrXII': 1078177,
            'chrXIII': 924431,
            'chrXIV': 784333,
            'chrXV': 1091291,
            'chrXVI': 948066
        }
        with open(wigfile) as fin:
            ga.add_from_wiggle(fin,"+")
            for chrom, length in chrdict.items():
                seg = GenomicSegment(chrom,0,length,"+")
                expected = ga[seg]
                found = self.bw[seg]
                yield self.check_vals_against_wig, expected, found
        
    def test_get_whole_chrom(self):
        assert False
        

    def test_iter(self):
        assert False