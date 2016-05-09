import os
import numpy

from nose.tools import assert_less_equal, assert_raises, assert_dict_equal,\
                       assert_true, assert_equal, assert_almost_equal

from pkg_resources import resource_filename
from plastid.readers.bigwig import BigWigReader
from plastid.readers.wiggle import WiggleReader
from plastid.genomics.roitools import GenomicSegment
from plastid.genomics.genome_array import GenomeArray
from plastid.util.services.decorators import skip_if_abstract


TOL=1e-5 #tolerance is high because bigwig files are approximate

base_path = resource_filename("plastid","test/data")
wigfile = os.path.join(base_path,"command_line","gen_reads_center_12_fw.wig")
bigwigfile = os.path.join(base_path,"command_line","gen_reads_center_12_fw.bw")


class WildCard(float):
    def __eq__(self,other):
        return True
    
    def __neq__(self,other):
        return False
    
    def __repr__(self):
        return 'any float'

    def __str__(self):
        return 'any float'

wildcard = WildCard()

class AbstractTestBBIFile():

    @classmethod
    def setUpClass(cls):
        cls.bw = BigWigReader(bigwigfile,fill=0) 

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
            
    @skip_if_abstract            
    def test_sum(self):
        assert False

#     @skip_if_abstract            
#     def test_summary_info(self):
#         assert False

class TestBigWigReader(AbstractTestBBIFile):
    
    @classmethod
    def setUpClass(cls):
        cls.bw = BigWigReader(bigwigfile,fill=0)
        cls.reader_class = BigWigReader 
        cls.chrdict = {
            'chrI'    : 230218,
            'chrII'   : 813184,
            'chrIII'  : 316620,
            'chrIV'   : 1531933,
            'chrV'    : 576874,
            'chrVII'  : 1090940,
            'chrVIII' : 562643,
            'chrIX'   : 439888,
            'chrX'    : 745751,
            'chrXI'   : 666816,
            'chrXII'  : 1078177,
            'chrXIII' : 924431,
            'chrXIV'  : 784333,
            'chrXV'   : 1091291,
            'chrXVI'  : 948066
        }
    
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
        with open(wigfile) as fin:
            ga.add_from_wiggle(fin,"+")
            for chrom, length in self.chrdict.items():
                seg = GenomicSegment(chrom,0,length,"+")
                expected = ga[seg]
                found = self.bw[seg]
                yield self.check_vals_against_wig, expected, found
    
    # NOTE: this test relies on WiggleReader being correct
    def check_random_windows_against_wig(self,strand):
        chrdict = self.chrdict
        chroms = list(self.chrdict)
        chridx = numpy.random.randint(0,high=len(chroms),size=50)
        ga = GenomeArray()
        
        i = 0
        
        with open(wigfile) as fin:
            ga.add_from_wiggle(fin,strand)
            while i < 50:
                chrom = chroms[chridx[i]]
                maxlength = chrdict[chrom]
                start = numpy.random.randint(0,high=maxlength-2000)
                end = numpy.random.randint(start+10000,high=start+20000)
                
                # make sure we don't go off chrom
                while end > maxlength:
                    end = numpy.random.randint(start+100,high=start+10000)
                    
                seg = GenomicSegment(chrom,start,end,strand)
                expected = ga[seg]
                # make sure segment has counts in it
                if expected.sum() > 0:
                    i += 1
                    found = self.bw[seg]
                    yield self.check_vals_against_wig, expected, found

    def test_random_windows_against_wig_fw(self):
        self.check_random_windows_against_wig("+")

    def test_random_windows_against_wig_rc(self):
        self.check_random_windows_against_wig("-")
    
    def test_get_chromosome_counts_zero_fill(self):
        ga = GenomeArray()
        with open(wigfile) as fin:
            ga.add_from_wiggle(fin,"+")
            for chrom, length in self.chrdict.items():
                seg = GenomicSegment(chrom,0,length,"+")
                expected = ga[seg]
                found = self.bw.get_chromosome_counts(chrom)
                yield self.check_vals_against_wig, expected, found

#     def test_get_chromosome_nan_fill(self):
#         bw = BigWigReader(bigwigfile,fill=numpy.nan)
#         ga = GenomeArray()
#         with open(wigfile) as fin:
#             ga.add_from_wiggle(fin,"+")
#             for chrom, length in self.chrdict.items():
#                 seg = GenomicSegment(chrom,0,length,"+")
#                 expected = ga[seg]
#                 found = bw.get_chromosome_counts(chrom)
#                 checkmask = ~numpy.isnan(found)
#                 
#                 # first to make sure no non-zero positions were masked
#                 # this doesn't guarantee we're ok, but will tell us if a big
#                 # problem is there
#                 yield assert_almost_equal, expected[checkmask].sum(),expected.sum()
#                 
#                 # there are no zeros in the input wiggle file; so all wiggles should
#                 # be nan
#                 yield assert_true, (~checkmask == (expected == 0)).all()
#                 
#                 # now check all non-nan positions
#                 yield self.check_vals_against_wig, expected[checkmask], found[checkmask]


    def test_fill_val_absent_chrom(self):
        filldef = BigWigReader(bigwigfile)
        fillnan = BigWigReader(bigwigfile,fill=numpy.nan)
        fill0   = BigWigReader(bigwigfile,fill=0)
        fill10  = BigWigReader(bigwigfile,fill=10)
        
        # chrVI is not in dataset; this should be an empty array
        seg = GenomicSegment("chrVI",5,1000,"+")
        
        assert_equal(len(filldef[seg]),len(seg),
                     "fetched wrong size")
        
#         assert_true(numpy.isnan(filldef[seg]).all(),
#                     "default not nan")
#         
#         assert_true(numpy.isnan(fillnan[seg]).all(),
#                     "nanfill didn't work")
        
        assert_true((fill0[seg] == 0).all(),
                    "0-fill didn't work")
        
#         assert_true((fill10[seg] == 10).all(),
#                     "10-fill didn't work")
        
    def test_fill_val_present_chrom(self):
        filldef = BigWigReader(bigwigfile)
        fillnan = BigWigReader(bigwigfile,fill=numpy.nan)
        fill0   = BigWigReader(bigwigfile,fill=0)
        fill10  = BigWigReader(bigwigfile,fill=10)
        
        # empty region
        seg = GenomicSegment("chrIV",5,10,"+")
        
        assert_equal(len(filldef[seg]),len(seg),
                     "fetched wrong size")
        
#         assert_true(numpy.isnan(filldef[seg]).all(),
#                     "default not nan")
#         
#         assert_true(numpy.isnan(fillnan[seg]).all(),
#                     "nanfill didn't work")
        
        assert_true((fill0[seg] == 0).all(),
                    "0-fill didn't work")
        
#         assert_true((fill10[seg] == 10).all(),
#                     "10-fill didn't work")
    
    def test_sum(self):
        bw = BigWigReader(os.path.join(base_path,"mini","wig","bw_fiveprime_15_fw.bw"))
        assert_equal(bw.sum(),4000)
    
    def test_iter(self):
        wig = WiggleReader(open(wigfile))
        bw  = BigWigReader(bigwigfile)
        
        for found in bw:
            expected = next(wig)
            fchrom = found[0]
            echrom = expected[0]
            assert_equal(fchrom,echrom,"Chromosome mismatch. Expected '%s', found '%s'." % (fchrom,echrom))

            fstart = found[1]
            estart = expected[1]
            assert_equal(fstart,estart,"Start mismatch. Expected '%s', found '%s'." % (fstart,estart))

            fend = found[2]
            eend = expected[2]
            assert_equal(fend,eend,"End mismatch. Expected '%s', found '%s'." % (fend,eend))

            fval = found[3]
            eval_ = expected[3]
            diff = abs(fval-eval_)
            assert_true(diff < TOL,"Difference %s exceeds tolerance '%s'. Expected '%s', found '%s'." % (diff,TOL,fval,eval_))

# Disabled until we decide what to do with summarize()            
#     def test_summarize(self):
#         bw  = BigWigReader(bigwigfile)
#         
#         #ga = GenomeArray(bw.chroms)
#         #ga.add_from_wiggle(open(wigfile),"+")
#         
#         chrom = "chrI"
#         maxlen = bw.chroms[chrom]
#         winstarts = numpy.random.randint(0,maxlen-20000,size=10)
#         winends   = winstarts + numpy.random.randint(500,40000,size=10)
#         winends[winends > maxlen] = maxlen
#         
#         numtests = 10
#         i = 0
#         while i < numtests:
#             s = numpy.random.randint(0,high=maxlen)
#             e = min(maxlen,numpy.random.randint(s+500,s+40000))
#             seg = GenomicSegment(chrom,s,e,"+")
#             #arr = ga[seg]
#             arr = bw[seg]
# 
#             labels = ["mean","max","min","cov","stdev"]
#             expected = [arr.mean(),arr.max(),arr.min(),wildcard,arr.std()]
#             
#             # change nans to 0
#             expected = [0 if numpy.isnan(X) else X for X in expected]
#             
#             print(expected)
#             found = bw.summarize(seg)
#             print(found)
#             print("---------------")
#             for label, exval, fval in zip(labels,expected,found):
#                 msg =  "test_summarize failed for stat '%s'. Expected %s, got %s (diff: %s)." % (label,exval,fval,abs(exval-fval))
#                 assert_almost_equal(exval,fval,msg=msg,delta=5)
#                 
#             i += 1
#         
#         # retval for summarize:  (mean,max_,min_,cov,stdev)
    
        
