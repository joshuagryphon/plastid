#!/usr/bin/env python
"""Tests and validates classes from :py:mod:`plastid.genomics.genome_array`,
these being |GenomeArray|, |SparseGenomeArray| and |BAMGenomeArray|,
using test data found in plastid.test.data. 

This module additionally contains utilites to generate other test datasets.
To do, please see the documentation for :py:func:`create_dataset`. Note,
this requires as external dependencies, Bowtie, Tophat, and Samtools.

"""
import copy
import tempfile
import os
import subprocess
import functools
import re
import unittest
import warnings

import pysam
import numpy.random

from pkg_resources import resource_filename, cleanup_resources
from nose.plugins.attrib import attr

from Bio import SeqIO


import plastid.util.services.exceptions

from plastid.readers.bed import BED_Reader
from plastid.genomics.genome_array import GenomeArray,\
                                       SparseGenomeArray,\
                                       BigWigGenomeArray,\
                                       BAMGenomeArray,\
                                       ThreePrimeMapFactory,\
                                       FivePrimeMapFactory,\
                                       CenterMapFactory,\
                                       five_prime_map,\
                                       three_prime_map,\
                                       center_map
from plastid.genomics.roitools import GenomicSegment, SegmentChain
from plastid.genomics.seqtools import random_seq
from plastid.util.io.openers import NullWriter
from plastid.util.services.decorators import skip_if_abstract
from plastid.util.services.mini2to3 import cStringIO


#===============================================================================
# INDEX: annotations/data used in unit tests and in generation of test datasets
#===============================================================================

# parameters to flesh out unit tests
# these are used by AbstractGenomeArrayHelper.set_class_parameters below
_GENOME_ARRAY_PARAMS = {
                        "test_class"    : GenomeArray,
                        "empty_regions" : ["splice","introns"],
                        "native_format" : "bowtie",
                        
                        }

_SPARSE_GENOME_ARRAY_PARAMS = {
                        "test_class"    : SparseGenomeArray,
                        "empty_regions" : ["splice","introns"],
                        "native_format" : "bowtie",
                        
                        }

_BIGWIG_GENOME_ARRAY_PARAMS = {
                        "test_class"    : BigWigGenomeArray,
                        "empty_regions" : ["introns"],
                        "native_format" : "bigwig",
                        
                        }

_BAM_GENOME_ARRAY_PARAMS = {
                        "test_class"    : BAMGenomeArray,
                        "empty_regions" : ["introns"],
                        "native_format" : "BAM",
                        
                        }

# descriptions of mapping configurations that we will use in test datasets
_SAMPLE_PAT = re.compile(r"(center|fiveprime|threeprime)_([0-9]+)")
_SAMPLE_BASES = ['center_0',
                 'center_12',
                 'fiveprime_0',
                 'fiveprime_15',
                 'threeprime_0',
                 'threeprime_15',
                 ]

_GA_MAP_FUNCTIONS = { "fiveprime"  : five_prime_map,
                      "threeprime" : three_prime_map,
                      "center"     : center_map,
                    }

_BAM_MAP_RULES = { "fiveprime_0"   : FivePrimeMapFactory(),
                   "threeprime_0"  : ThreePrimeMapFactory(),
                   "center_0"      : CenterMapFactory(),
                   "fiveprime_15"  : FivePrimeMapFactory(15),
                   "threeprime_15" : ThreePrimeMapFactory(15),
                   "center_12"     : CenterMapFactory(12),
                  }

# constants/names of files in test datasets, to use in test cases 
# or to generate using methods below
#
# all filenames are relative to a base_folder that is passed to individual functions
# these files will be created by create_test_dataset(), below
_TEST_FILES = { "variable_step_fw"    : os.path.join("wig","variable_step_fw.wig"),
                "variable_step_rc"    : os.path.join("wig","variable_step_rc.wig"),
                "bedgraph_fw"  : os.path.join("wig","bedgraph_fw.wig"),
                "bedgraph_rc"  : os.path.join("wig","bedgraph_rc.wig"),
                "juncs"        : os.path.join("ebwt","chrA.juncs"),
                "bowtie_index" : os.path.join("ebwt","chrA"),
                "bed"          : os.path.join("bed","chrA.bed"),
                "reads"        : os.path.join("fasta","chrA_reads.fa"),
                "bowtie"       : os.path.join("align","chrA_unspliced.bowtie"),
                "bam"          : os.path.join("align","chrA_tophat.bam"),
                "genome"       : os.path.join("fasta","chrA.fa"),
                }

# annotation data
TEST_CHR_BED="""chrA    100    1100    unique_plus    0    +    -1    -1    0,0,0    1    1000,    0,
chrA    100    1100    unique_minus    0    -    -1    -1    0,0,0    1    1000,    0,
chrA    1200    2250    entire_repeat_region_plus    0    +    -1    -1    0,0,0    1    1050,    0,
chrA    1200    2250    entire_repeat_region_minus    0    -    -1    -1    0,0,0    1    1050,    0,
chrA    1200    1700    repeat_1_plus    0    +    -1    -1    0,0,0    1    500,    0,
chrA    1200    1700    repeat_1_minus    0    -    -1    -1    0,0,0    1    500,    0,
chrA    1750    2250    repeat_2_plus    0    +    -1    -1    0,0,0    1    500,    0,
chrA    1750    2250    repeat_2_minus    0    -    -1    -1    0,0,0    1    500,    0,
chrA    2350    2475    splice_plus    100    +    -1    -1    0,0,0    2    25,25,    0,100,
chrA    2350    2475    splice_minus    100    -    -1    -1    0,0,0    2    25,25,    0,100,
chrA    2375    2450    intron_plus    0    +    -1    -1    0,0,0    1    75,    0,
chrA    2375    2450    intron_minus    0    -    -1    -1    0,0,0    1    75,    0,""".replace("    ","\t")

TEST_CHR_JUNCS="""chrA    2374    2450    +
chrA    2374    2450    -""".replace("    ","\t")

# miscellaneous constants
STRAND_KEYS = { "+" : "fw", "-" : "rc" }
DEFAULT_READS_PER_REGION = 1000
DEFAULT_READ_LENGTH = 30


#===============================================================================
# INDEX: Helper functions for unit tests and test dataset creation methods below
#===============================================================================

def tearDownModule():
    """Remove test dataset files after unit tests are complete"""
    cleanup_resources()

def fetch_regions():
    """Parses test regions of interest for synthetic genomes
    
    Returns
    -------
    list<SegmentChain>
    """
    return list(BED_Reader(cStringIO.StringIO(TEST_CHR_BED),return_type=SegmentChain))

def _read_count_vectors(base_folder):
    """Read count vectors from a synthetic datasets
    generated by :py:method:`create_test_dataset`

    Parameters
    ----------
    base_folder : str
        path to base folder passed to :py:method:`create_test_dataset`
    
    Returns
    -------
    dict : dict of numpy.ndarrays of count data
    """
    dtmp = {}
    for k in _SAMPLE_BASES:
        for strand_key in ("fw","rc"):
            dtmp["%s_%s" % (k,strand_key) ] = numpy.loadtxt(os.path.join(base_folder,"count_vectors","%s_%s.txt" % (k,strand_key)))
    return dtmp

def _read_bowtie_files_to_genome_arrays(base_folder,test_class=GenomeArray):
    """Construct |GenomeArray| s from bowtie files

    Parameters
    ----------
    base_folder : str
        path to base folder passed to create_test_dataset()

    test_class : class
        Subclass of |MutableGenomeArray| (e.g. |GenomeArray| or |SparseGenomeArray| to test)

    Returns
    -------
    dict : dict of |GenomeArray| s of mapped read alignments from bowtie
    """    
    gnds = {}
    for k in _SAMPLE_BASES:
        mapping, offset = _SAMPLE_PAT.search(k).groups()
        trans_key = "nibble" if mapping == "center" else "offset"
        trans_args = { trans_key : int(offset) }
        gnds[k] = test_class()
        gnds[k].add_from_bowtie(open(os.path.join(base_folder,_TEST_FILES["bowtie"])),_GA_MAP_FUNCTIONS[k.split("_")[0]],**trans_args)
    
    return gnds

def _get_ivc_numpy_counts(ivc,count_vec):
    """Fetches appropriately-spliced counts at each position in an ROI from a numpy array
    
    Parameters
    ----------
    ivc : |SegmentChain|
        SegmentChain describing region of interest
    
    count_vec : numpy.ndarray
        numpy.ndarray, in coordinates matching those of ivc
    
    Returns
    -------
    numpy.ndarray : numpy.ndarray of counts each each position in ivc
    """
    counts = []
    for iv in ivc:
        counts.extend(count_vec[iv.start:iv.end])
    
    if ivc.spanning_segment.strand == "-":
        counts = counts[::-1]
    
    return numpy.array(counts)


#===============================================================================
# INDEX: unittest suites
#===============================================================================
        
class AbstractGenomeArrayHelper(unittest.TestCase):
    """Abstract base class for various types of |AbstractGenomeArray| test cases"""

    set_up = False
    
    @staticmethod
    def set_class_parameters(cls,params,test_folder=resource_filename("plastid","test/data/mini"),tol=1e-8):
        """Set class parameters on the creation of the first instance.
        This is a bit of a hack because we need to set class parameters.
        We can't do this in a ``setUpClass`` method, because ``setUpClass`` only 
        accepts a single parameter (the class). We don't want to do this in
        ``__init__`` either, because unittest calls ``__init__`` once per test run,
        and these operations are expensive. So, instead we define this method,
        and call it from ``__init__`` if and only if ``cls.set_up == False``
        
        Parameters
        ----------
        cls : class
            class that is a subclass of :py:class:`unittest.TestCase`,
            to which parameters will be appended
        
        params : dict
            Parameters specific to the set-up of test suites for specific
            types of GenomeArrays
        
        test_folder : str or :py:class:`Resource`
            Real or virtual location of folder of test data

        tol : float
            Tolerance for numerical differences between expected and observed
            values in the various tests        
        """
        cls.test_folder = test_folder
        cls.tol = tol
        
        cls.test_class = params["test_class"]
        cls.native_format = params["native_format"]
        cls.count_vecs = _read_count_vectors(cls.test_folder)
        
        cls.regions = fetch_regions()
        cls.region_classes = {
                                "unique"  : [X for X in cls.regions if "unique" in X.get_name()],
                                "repeat"  : [X for X in cls.regions if "entire" not in X.get_name() and "repeat" in X.get_name()],
                                "introns" : [X for X in cls.regions if "intron" in X.get_name()],
                                "splice"  : [X for X in cls.regions if "splice" in X.get_name()],
                                "entire"  : [X for X in cls.regions if "entire" in X.get_name()],
                               }
        
        cls.region_classes["empty"] = []
        cls.empty_names = []
        for k in params["empty_regions"]:
            my_regions = cls.region_classes[k]
            cls.region_classes["empty"].extend(my_regions)
            cls.empty_names.extend([X.get_name() for X in my_regions])
        
        cls.expected_unnorm_sum = 0
        #for region in set(cls.regions) - set(cls.region_classes["empty"]) - set(cls.region_classes["entire"]):
        read_regions = [X for X in cls.regions if all([X.get_name() not in cls.empty_names,"entire" not in X.get_name()])]
        for region in read_regions:
            vec_key = "fw" if region.strand == "+" else "rc"
            cls.expected_unnorm_sum += _get_ivc_numpy_counts(region,cls.count_vecs["fiveprime_0_%s" % vec_key]).sum()
        
        cls.set_up = True
        
    def __init__(self,
                 methodName='runTest',
                 params={},
                 test_folder=resource_filename("plastid","test/data/mini"),
                 tol=1e-8):
        """Initialize test case to run a single method. 
        We override this method to make sure expensive operations are only run when
        the first instance is made, and then stored in class attributes

        Parameters
        ----------
        methodName : str
            Name of method being run. Required by :py:class:`unittest.TestCase`

        params : dict
            Parameters specific to the set-up of test suites for specific
            |AbstractGenomeArray| subclasses. Don't change these
            
        test_folder : str or :py:class:`Resource`
            Real or virtual location of folder of test data

        tol : float
            Tolerance for numerical differences between expected and observed
            values in the various tests
        """
        unittest.TestCase.__init__(self,methodName=methodName)
        # only do setup if __init__ is called by a subclass
        if "Abstract" not in self.__class__.__name__:
            if self.__class__.set_up == False:
                AbstractGenomeArrayHelper.set_class_parameters(self.__class__,
                                                               params=params,
                                                               test_folder=test_folder,
                                                               tol=tol)

    @skip_if_abstract
    def test_chroms(self):
        for v in self.gnds.values():
            self.assertEqual(set(v.chroms()),set(["chrA"]))
        
    @skip_if_abstract
    def test_strands(self):
        possible_strands = set(["+","-","."])
        for v in self.gnds.values():
            self.assertGreaterEqual(len(set(v.strands()) & possible_strands),0)
        
    @skip_if_abstract
    def test_test_class(self):
        # Assure all genome arrays tested are of correct subclass
        for k,v in self.gnds.items():
            self.assertTrue(isinstance(v,self.test_class),
                            "Test %s: instance is of wrong class (expected: %s, found %s)" % (k,
                            self.test_class.__name__,
                            v.__class__.__name__))
 
    @skip_if_abstract
    def test_native_import_positionwise_equality_unique_regions(self):
        for k in _SAMPLE_BASES:
            for region in self.region_classes["unique"]:
                strand_key = STRAND_KEYS[region.spanning_segment.strand]
                gnd_counts   = numpy.array(region.get_counts(self.gnds[k]))
                self.assertGreater(gnd_counts.sum(),0,"Region is empty in sample %s" % k)
                known_counts = _get_ivc_numpy_counts(region,self.count_vecs["%s_%s" % (k,strand_key)])
                
                max_err = max(abs(gnd_counts - known_counts))
                msg1 = "Positionwise count difference (%s) exceeded tolerance (%s) for '%s' file import for sample test '%s'" % (self.tol,max_err,self.native_format,k)
                self.assertLessEqual(max_err,self.tol,msg1)
                
                sum_diff = abs(known_counts.sum() - gnd_counts.sum())
                msg2 = "Error in difference of total counts (%s) exceeded tolerance (%s) for '%s' import for sample test %s" % (sum_diff,self.tol,self.native_format,k)
                self.assertLessEqual(sum_diff,self.tol,msg2)
        
    @skip_if_abstract
    def test_native_import_positionwise_equality_repeat_regions(self):
        # test sums of position-wise vectors for repeat regions
        plus_repeat  = [X for X in self.region_classes["repeat"] if X.spanning_segment.strand == "+"]
        minus_repeat = [X for X in self.region_classes["repeat"] if X.spanning_segment.strand == "-"]

        
        lengths = set([X.length for X in plus_repeat + minus_repeat])
        self.assertEqual(len(lengths),1)

        for k in _SAMPLE_BASES:
            plus_vec  = numpy.zeros(plus_repeat[0].length)
            minus_vec = numpy.zeros(plus_repeat[0].length)
            
            known_plus_vec  = numpy.zeros(plus_repeat[0].length)
            known_minus_vec = numpy.zeros(plus_repeat[0].length)
            
            for region in plus_repeat:
                plus_vec += region.get_counts(self.gnds[k])
                known_plus_vec += _get_ivc_numpy_counts(region,
                                                        self.count_vecs["%s_%s" % (k,"fw")])
    
            for region in minus_repeat:
                minus_vec += region.get_counts(self.gnds[k])
                known_minus_vec += _get_ivc_numpy_counts(region,
                                                         self.count_vecs["%s_%s" % (k,"rc")])
            
            self.assertGreater(plus_vec.sum(),0)
            self.assertGreater(minus_vec.sum(),0)
    
            self.assertTrue((abs(known_plus_vec - plus_vec)<=self.tol).all(),
                            "Positionwise count difference exceeded tolerance %s for %s import on sample test %s on plus strand" % (self.tol,self.native_format,k))
            self.assertTrue((abs(known_minus_vec - minus_vec)<=self.tol).all(),
                            "Positionwise count difference exceeded tolerance %s for %s import for sample test %s on minus strand" % (self.tol,self.native_format,k))

    @skip_if_abstract
    def test_native_import_empty_regions(self):
        # test regions that should be empty (e.g. introns and splicing)
        for k in _SAMPLE_BASES:
            for region in self.region_classes["empty"]:
                self.assertEqual(sum(region.get_counts(self.gnds[k])),0,
                                 "Found counts in region that should be empty for sample test %s" % k)

    @skip_if_abstract
    def variablestep_and_bedgraph_export_helper(self,wiggle_type,export_function,input_class=None,**kwargs):
        """Helper function to evaluate tests on variable step wiggle or BED export

        Parameters
        ----------
        wiggle_type : str
            Type of wiggle file. "variable_step" or "bedgraph"

        export_function : function
            unbound method defining export type (e.g. GenomeArray.to_variable_step, BAMGenomeArray.to_bedgraph)

        input_class : subclass of |MutableAbstractGenomeArray| or None
            Class into which exported wiggle or bedgraph files will be read. If None, defaults to self.test_class

        kwargs : keyword arguments
        """
        if input_class is None:
            input_class = self.test_class

        for k,v in self.gnds.items():
            fw_out = tempfile.NamedTemporaryFile(mode="w",delete=False)
            rc_out = tempfile.NamedTemporaryFile(mode="w",delete=False)
            
            export_function(v,fw_out,"test","+",**kwargs)
            export_function(v,rc_out,"test","-",**kwargs)
    
            fw_out.close()
            rc_out.close()
            
            new_gnd = input_class()
            new_gnd.add_from_wiggle(open(fw_out.name),"+")
            new_gnd.add_from_wiggle(open(rc_out.name),"-")

            self.assertGreater(v.lengths()["chrA"],0)
            self.assertGreater(new_gnd.lengths()["chrA"],0)
            
            ivplus  = GenomicSegment("chrA",0,v.lengths()["chrA"],"+")
            ivminus = GenomicSegment("chrA",0,v.lengths()["chrA"],"-")

            # test equality of what was exported with current state of GenomeArray
            self.assertTrue(abs(new_gnd[ivplus] - v[ivplus] <= self.tol).all(),
                            "%s wiggle output on plus strand failed positionwise tolerance %s for test %s" % (wiggle_type,self.tol,k))
            self.assertGreater(new_gnd[ivplus].sum(),0,
                               "No counts found for %s reimport test %s" % (wiggle_type,k))

            self.assertTrue(abs(new_gnd[ivminus] - v[ivminus] <= self.tol).all(),
                            "%s wiggle output on minus strand failed positionwise tolerance %s for test %s" % (wiggle_type,self.tol,k))
            self.assertGreater(new_gnd[ivminus].sum(),0,
                               "No counts found for %s reimport test %s" % (wiggle_type,k))
            
            # ground-truth test against numpy arrays for unique regions
            for region in self.region_classes["unique"]:
                strand_key = STRAND_KEYS[region.spanning_segment.strand]
                gnd_counts   = numpy.array(region.get_counts(new_gnd))
                self.assertGreater(gnd_counts.sum(),0,"Reimported region is empty in sample %s" % k)
                known_counts = _get_ivc_numpy_counts(region,self.count_vecs["%s_%s" % (k,strand_key)])
                max_err = max(abs(gnd_counts - known_counts))
                self.assertLessEqual(max_err,self.tol,
                                "Positionwise count difference (%s) exceeded tolerance (%s) for %s reimport after export from class '%s' for sample '%s'" % (self.tol,max_err,self.native_format,self.test_class,k))

            os.remove(fw_out.name)
            os.remove(rc_out.name)

    @skip_if_abstract
    def test_unnormalized_sum(self):
        for k,v in self.gnds.items():
            v.set_normalize(False)
            found = v.sum()
            expected = self.expected_unnorm_sum
            err = abs(found - expected)
            err_msg = "Observed error (%s) in unnormalized sum (observed %s; expected %s) greater than tolerance (%s) for sample '%s'" % (err,found,expected,self.tol,k)
            self.assertLessEqual(err,self.tol,err_msg)
    
    @skip_if_abstract
    def test_normalize_not_change_sum(self):
        for k,v in self.gnds.items():
            v.set_normalize(True)
            found_sum = v.sum()
            err_msg = "Normalize flag changed sum to %s from %s for sample '%s'" % (found_sum,self.expected_unnorm_sum,k)
            err = abs(found_sum - self.expected_unnorm_sum)
            self.assertLessEqual(err,self.tol,err_msg)
            v.set_normalize(False)

    @skip_if_abstract
    def test_set_and_reset_sum(self):
        expected_unnorm_sum2 = 50000
        for k,v in self.gnds.items():
            v.set_sum(expected_unnorm_sum2)
            v.set_normalize(False)
            
            found = v.sum()
            err = abs(found - expected_unnorm_sum2)
            err_msg = "Observed error (%s) in sample '%s' set unnormalized sum (observed %s; expected %s) greater than tolerance %s" % (err,k,found,expected_unnorm_sum2,self.tol)
            self.assertLessEqual(err,self.tol,err_msg)

            v.set_normalize(True)
            found = v.sum()
            err = abs(found - expected_unnorm_sum2)
            err_msg = "Observed error (%s) in sample '%s' set normalized sum (observed %s; expected %s) greater than tolerance %s" % (err,k,found,expected_unnorm_sum2,self.tol)
            self.assertLessEqual(err,self.tol,err_msg)

            v.set_normalize(False)
            v.reset_sum()
            found = v.sum()
            err = abs(found - self.expected_unnorm_sum)
            err_msg = "Observed error (%s) in sample '%s' reset sum (observed %s; expected %s) greater than tolerance %s" % (err,k,found,self.expected_unnorm_sum,self.tol)
            self.assertLessEqual(err,self.tol,err_msg)

    @skip_if_abstract
    def test_regionwise_normalize_and_sum(self):
        expected_unnorm_sum2 = 50000
        for k, v in self.gnds.items():
            v.set_normalize(False)
            v.reset_sum()
            
            # add an order of magnitude to account for summing
            tol = self.tol*10

            # exclude repeat regions, because those will align differently than they were generated
            # remove "empty" also, because this will include spliced regions for some tests, 
            # as necessary
            nonrepeat_nonempty = [X for X in self.regions if all(["repeat" not in X.get_name(),X.get_name() not in self.empty_names])]
            for region in nonrepeat_nonempty: #set(self.regions) - set(self.region_classes["repeat"]) - set(self.region_classes["empty"]):
                # Make sure baseline number is ok
                found_region_sum = sum(region.get_counts(v))
                expected_region_unnorm = _get_ivc_numpy_counts(region,self.count_vecs["%s_%s" % (k,STRAND_KEYS[region.strand])]).sum()
                err = abs(found_region_sum-expected_region_unnorm)
                self.assertLessEqual(err,tol,
                                     "Found unnormalized region sum %s different from expected %s more than error %s" % (found_region_sum,expected_region_unnorm,tol))

                # Test normalize
                v.set_normalize(True)
                expected_region_norm  = float(expected_region_unnorm) / v.sum() * 10**6
                found_region_sum = sum(region.get_counts(v))
                err = abs(found_region_sum-expected_region_norm)
                self.assertLessEqual(err,tol,
                                     "Found normalized region sum (%s) different from expected (%s) more than error (observed %s; tolerance %s) for sample '%s'" % (found_region_sum,expected_region_norm,err,tol,k))

                # Test reversibility
                v.set_normalize(False)
                found_region_sum = sum(region.get_counts(v))
                err = abs(found_region_sum-expected_region_unnorm)
                self.assertLessEqual(err,tol,
                                     "Found re-unnormalized region sum %s different from expected %s more than error %s" % (found_region_sum,expected_region_unnorm,tol))

                # Set sum, no normalization
                v.set_sum(expected_unnorm_sum2)
                found_region_sum = sum(region.get_counts(v))
                err = abs(found_region_sum-expected_region_unnorm)
                self.assertLessEqual(err,tol,
                                     "Found post-global-sum-set unnormalized region sum %s different from expected %s more than error %s" % (found_region_sum,expected_region_unnorm,tol))

                # Add normalization on top of set sum
                v.set_normalize(True)
                expected_region_norm2 = float(expected_region_unnorm) / expected_unnorm_sum2 * 10**6
                found_region_sum = sum(region.get_counts(v))
                err = abs(found_region_sum-expected_region_norm2)
                self.assertLessEqual(err,tol,
                                     "Found post-global-sum-set normalized region sum %s different from expected %s more than error %s" % (found_region_sum,expected_region_norm2,tol))

                # Reset sum, keep normalization
                v.reset_sum()
                found_region_sum = sum(region.get_counts(v))
                err = abs(found_region_sum-expected_region_norm)
                self.assertLessEqual(err,tol,
                                     "Found post-reset normalized region sum %s different from expected %s more than error %s" % (found_region_sum,expected_region_norm,tol))

                # Revert all
                v.set_normalize(False)
                found_region_sum = sum(region.get_counts(v))
                err = abs(found_region_sum-expected_region_unnorm)
                self.assertLessEqual(err,tol,
                                     "Found unnormalized region sum %s different from expected %s more than error %s" % (found_region_sum,expected_region_unnorm,tol))

    @skip_if_abstract
    def test_get_genomicsegment_roi_order_false(self):
        k = _SAMPLE_BASES[0]
        for region in self.region_classes["unique"]:
            seg = region.spanning_segment
            strand_key = STRAND_KEYS[region.spanning_segment.strand]
            #gnd_counts   = self.gnds[k].__getitem__(seg,roi_order=False)
            gnd_counts   = self.gnds[k].get(seg,roi_order=False)
            known_counts = self.count_vecs["%s_%s" % (k,strand_key)][seg.start:seg.end]

            max_err = max(abs(gnd_counts - known_counts))
            self.assertLessEqual(max_err,self.tol,
                            "Positionwise count difference '%s' exceeded tolerance '%s' for %s __getitem__ with roi_order==False for sample test %s" % (self.tol,max_err,self.native_format,k))


    @skip_if_abstract
    def test_get_genomicsegment_roi_order_true(self):
        k = _SAMPLE_BASES[0]
        for region in self.region_classes["unique"]:
            seg = region.spanning_segment
            strand_key = STRAND_KEYS[region.spanning_segment.strand]
            gnd_counts   = self.gnds[k].get(seg,roi_order=True)
            known_counts = self.count_vecs["%s_%s" % (k,strand_key)][seg.start:seg.end]
            if seg.strand == "-":
                known_counts = known_counts[::-1]

            max_err = max(abs(gnd_counts - known_counts))
            self.assertLessEqual(max_err,self.tol,
                            "Positionwise count difference '%s' exceeded tolerance '%s' for %s __getitem__ with roi_order==True for sample test %s" % (self.tol,max_err,self.native_format,k))
    @skip_if_abstract
    def test_getitem_genomicsegment(self):
        k = _SAMPLE_BASES[0]
        for region in self.region_classes["unique"]:
            seg = region.spanning_segment
            strand_key = STRAND_KEYS[region.spanning_segment.strand]
            gnd_counts   = self.gnds[k].__getitem__(seg)
            known_counts = self.count_vecs["%s_%s" % (k,strand_key)][seg.start:seg.end]
            if seg.strand == "-":
                known_counts = known_counts[::-1]

            max_err = max(abs(gnd_counts - known_counts))
            self.assertLessEqual(max_err,self.tol,
                            "Positionwise count difference '%s' exceeded tolerance '%s' for %s __getitem__ with roi_order==True for sample test %s" % (self.tol,max_err,self.native_format,k))

    @skip_if_abstract
    def test_getitem_segmentchain(self):
        k = _SAMPLE_BASES[0]
        for region in self.region_classes["unique"]:
            strand_key = STRAND_KEYS[region.spanning_segment.strand]
            gnd_counts   = self.gnds[k][region] # test
            self.assertGreater(gnd_counts.sum(),0,"Region is empty in sample %s" % k)
            known_counts = _get_ivc_numpy_counts(region,self.count_vecs["%s_%s" % (k,strand_key)])
            max_err = max(abs(gnd_counts - known_counts))
            self.assertLessEqual(max_err,self.tol,
                            "Positionwise count difference '%s' exceeded tolerance '%s' for %s import for sample test %s" % (self.tol,max_err,self.native_format,k))



class AbstractExportableGenomeArrayHelper(AbstractGenomeArrayHelper):
    @skip_if_abstract
    def test_variablestep_export(self):
        self.variablestep_and_bedgraph_export_helper("variable_step",self.test_class.to_variable_step)

    @skip_if_abstract
    def test_bedgraph_export(self):
        self.variablestep_and_bedgraph_export_helper("bedgraph",self.test_class.to_bedgraph)


@attr(test="unit")
@attr(speed="slow")
class TestGenomeArray(AbstractExportableGenomeArrayHelper):
    """Test case for :py:class:`GenomeArray`"""

    set_up   = False
    has_gnds = False
    
    def __init__(self,
                 methodName='runTest',
                 params=_GENOME_ARRAY_PARAMS,
                 test_folder=resource_filename("plastid","test/data/mini"),
                 tol=1e-8):
        """Initialize test case to run a single method. 
        We override this method to make sure expensive operations are only run when
        the first instance is made, and then stored in class attributes

        Parameters
        ----------
        methodName : str
            Name of method being run. Required by :py:class:`unittest.TestCase`
        
        params : dict
            Parameters specific to the set-up of test suites for specific
            |AbstractGenomeArray| subclasses. Don't change these
            
        test_folder : str or Resource
            Real or virtual location of folder of test data

        tol : float
            Tolerance for numerical differences between expected and observed
            values in the various tests
        """
        AbstractGenomeArrayHelper.__init__(self,methodName=methodName,params=params,test_folder=test_folder,tol=tol)
        if self.__class__.has_gnds == False:
            self.__class__.gnds = _read_bowtie_files_to_genome_arrays(self.test_folder,self.test_class)
            self.__class__.has_gnds = True
            #TestGenomeArray.setUpClassOnlyOnce()

    def test_setitem_genomicsegment_scalar(self):
        ga = GenomeArray({"chrA" : 2000})
        segplus  = GenomicSegment("chrA",50,100,"+")
        segminus = GenomicSegment("chrA",50,100,"-")

        # scalar set
        ga[segplus]  = 52
        ga[segminus] = 342

        self.assertTrue((ga._chroms["chrA"]["+"][50:100]==52).all(),"%s failed scalar genomicsegment __setitem__ for plus strand.")
        self.assertTrue((ga._chroms["chrA"]["-"][50:100]==342).all(),"%s failed scalar genomicsegment __setitem__ for minus strand.")
        self.assertEqual(ga.sum(),52*len(segplus) + 342*len(segminus))

    def test_setitem_genomicsegment_vector(self):
        ga = GenomeArray({"chrA" : 2000})
        segplus  = GenomicSegment("chrA",50,100,"+")
        segminus = GenomicSegment("chrA",50,100,"-")

        # vector set
        r1 = numpy.random.randint(0,high=242,size=50)
        r2 = numpy.random.randint(0,high=242,size=50)
        ga[segplus]  = r1
        ga[segminus] = r2

        self.assertTrue((ga._chroms["chrA"]["+"][50:100]==r1).all(),"%s failed vector genomicsegment __setitem__ for plus strand.")
        self.assertTrue((ga._chroms["chrA"]["-"][50:100]==r2[::-1]).all(),"%s failed vector genomicsegment __setitem__ for minus strand.")

        self.assertEqual(ga.sum(),r1.sum()+r2.sum())

    def test_setitem_segmentchain_scalar(self):
        ga = GenomeArray({"chrA" : 2000})
        pluschain = SegmentChain(GenomicSegment("chrA",50,100,"+"),
                                 GenomicSegment("chrA",150,732,"+"),
                                 GenomicSegment("chrA",1800,2500,"+"))
        minuschain = SegmentChain(GenomicSegment("chrA",50,100,"-"),
                                  GenomicSegment("chrA",150,732,"-"),
                                  GenomicSegment("chrA",1800,2500,"-"))
        ga[pluschain]  = 31
        ga[minuschain] = 424
        for seg in pluschain:
            self.assertTrue((ga._chroms[seg.chrom][seg.strand][seg.start:seg.end]==31).all())

        for seg in minuschain:
            self.assertTrue((ga._chroms[seg.chrom][seg.strand][seg.start:seg.end]==424).all())
        
        self.assertEqual(ga.sum(),31*pluschain.length+424*minuschain.length)

    def test_setitem_segmentchain_vector(self):
        ga = GenomeArray({"chrA" : 2000})
        pluschain = SegmentChain(GenomicSegment("chrA",50,100,"+"),
                                 GenomicSegment("chrA",150,732,"+"),
                                 GenomicSegment("chrA",1800,2500,"+"))
        minuschain = SegmentChain(GenomicSegment("chrA",50,100,"-"),
                                  GenomicSegment("chrA",150,732,"-"),
                                  GenomicSegment("chrA",1800,2500,"-"))

        plusvec  = numpy.random.randint(0,high=250,size=pluschain.length)
        minusvec = numpy.random.randint(0,high=250,size=minuschain.length)

        ga[pluschain]  = plusvec
        ga[minuschain] = minusvec

        x = 0
        for seg in pluschain:
            subvec = ga._chroms["chrA"]["+"][seg.start:seg.end]
            self.assertTrue((subvec==plusvec[x:x+len(subvec)]).all())
            x += len(subvec)

        x = 0
        for seg in minuschain:
            subvec = ga._chroms["chrA"]["-"][seg.start:seg.end][::-1]
            self.assertTrue((subvec==minusvec[len(minusvec)-x-len(subvec):len(minusvec)-x]).all())
            x += len(subvec)

        self.assertEqual(ga.sum(),plusvec.sum()+minusvec.sum())

    def variablestep_and_bed_import_helper(self,wiggle_type):
        """Helper function to evaluate tests on variable step wiggle or BEDgraph import

        Parameters
        ----------
        wiggle_type : str
            Type of wiggle file. "variable_step" or "bedgraph"
        """
        gnd = self.test_class()
        gnd.add_from_wiggle(open(os.path.join(self.test_folder,_TEST_FILES["%s_%s" % (wiggle_type,"fw")])),"+")
        gnd.add_from_wiggle(open(os.path.join(self.test_folder,_TEST_FILES["%s_%s" % (wiggle_type,"rc")])),"-")

        # Make sure imported counts are nonzero
        self.assertGreater(gnd.sum(),0,"Import of %s yielded no counts!" % wiggle_type)
        chrA_len = gnd.lengths()["chrA"]
        
        for strand,trackstub,label in [("+","fw","plus"),("-","rc","minus")]:
            my_vec = self.count_vecs["fiveprime_0_%s" % trackstub]
            vec_len = len(my_vec)

            empty_iv = GenomicSegment("chrA",vec_len,chrA_len,strand)
            nonempty_iv = GenomicSegment("chrA",0,vec_len,strand)
            nonempty_vec = gnd.get(nonempty_iv,roi_order=False)

            # make sure count vector has requisite counts
            self.assertGreater(my_vec.sum(),0)

            # make sure sums are equal
            self.assertEquals(my_vec.sum(),nonempty_vec.sum())

            # make sure all regions after count vector are empty
            self.assertEquals(gnd[empty_iv].sum(),0,"Found counts in region that should be empty.")

            # make sure all positions in regions up to count vector are correct
            max_err = abs(my_vec - nonempty_vec).max()
            self.assertLessEqual(max_err,self.tol,
                                 "Positionwise count difference %s exceeded tolerance %s for wiggle import on %s strand" % (max_err,self.tol,label)
                                 )
        
    def test_bedgraph_import(self):
        self.variablestep_and_bed_import_helper("bedgraph")
    
    def test_variablestep_import(self):
        self.variablestep_and_bed_import_helper("variable_step")
   
    def test_genome_wide_scalar_plus_equals_sum(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        expected_length = 2*sum(chroms.values())
        
        gnd = self.test_class(chroms)
        self.assertEqual(len(gnd),expected_length)
        
        gnd += 5
        self.assertEqual(gnd.sum(),5*expected_length)
        
        gnd -= 5
        self.assertEqual(gnd.sum(),0)
    
    def test_genome_wide_scalar_plus_equals_not_change_length(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        expected_length = 2*sum(chroms.values())
        
        gnd = self.test_class(chroms)
        self.assertEqual(len(gnd),expected_length)
        
        gnd += 5
        self.assertEqual(len(gnd),expected_length)

    def test_genome_wide_scalar_plus_sum(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        expected_length = 2*sum(chroms.values())
        
        gnd = self.test_class(chroms)
        self.assertEqual(len(gnd),expected_length)
        
        # add scalar
        gnd2 = gnd + 5
        self.assertEqual(gnd2.sum(),5*expected_length)
        self.assertEqual(gnd.sum(),0)
        
        # add scalar to occupied gnd
        gnd3 = gnd2 + 1
        self.assertEqual(gnd3.sum(),6*expected_length)

    def test_genome_scalar_times_equals_sum(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        expected_length = 2*sum(chroms.values())
        
        gnd = self.test_class(chroms)
        self.assertEqual(len(gnd),expected_length)
        
        gnd *= 5
        self.assertEqual(gnd.sum(),0)
        
        gnd += 1
        gnd *= 5
        self.assertEqual(gnd.sum(),5*expected_length)
    
    def test_genome_wide_scalar_times_equals_not_change_length(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        expected_length = 2*sum(chroms.values())
        
        gnd = self.test_class(chroms)
        self.assertEqual(len(gnd),expected_length)
        
        gnd *= 5
        self.assertEqual(len(gnd),expected_length)

    def test_genome_wide_scalar_times(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        expected_length = 2*sum(chroms.values())
        
        gnd  = self.test_class(chroms)
        gnd += 1
        self.assertEqual(gnd.sum(),expected_length)
        
        gnd2 = gnd*2
        self.assertEqual(gnd2.sum(),2*expected_length)
        self.assertEqual(gnd.sum(),expected_length)
   
    def test_genome_wide_array_add_same_size(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        
        gnd1 = self.test_class(chroms)
        gnd2 = self.test_class(chroms)
        
        gnd1 += 1
        gnd2 += 3
        gnd3 = gnd1 + gnd2
        self.assertEquals(gnd3.sum(),gnd2.sum()+gnd1.sum())
        
        iv1plus  = GenomicSegment("chrA",0,1000,"+")
        iv1minus = GenomicSegment("chrA",0,1000,"+")
        
        self.assertTrue((gnd3[iv1plus]==gnd2[iv1plus]+gnd1[iv1plus]).all())
        self.assertTrue((gnd3[iv1minus]==gnd2[iv1minus]+gnd1[iv1minus]).all())

    def test_genome_wide_array_multiply_same_size(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        
        gnd1 = self.test_class(chroms)
        gnd2 = self.test_class(chroms)
        
        gnd1 += 2
        gnd2 += 3
        gnd3 = gnd1 * gnd2

        iv1plus  = GenomicSegment("chrA",0,1000,"+")
        iv1minus = GenomicSegment("chrA",0,1000,"+")
        
        self.assertTrue((gnd3[iv1plus]==gnd2[iv1plus]*gnd1[iv1plus]).all())
        self.assertTrue((gnd3[iv1minus]==gnd2[iv1minus]*gnd1[iv1minus]).all())
   
    def test_iadd_no_normalize(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        gnd = self.test_class(chroms)
        gnd.set_normalize(False)
        self.assertEqual(0,gnd.sum())
        
        iv1plus  = GenomicSegment("chrA",0,1000,"+")
        iv2plus  = GenomicSegment("chrA",0,500,"+")
        gnd[iv1plus] += 1
        self.assertEqual(1000,gnd.sum())

        gnd[iv1plus] += 1
        self.assertEqual(2000,gnd.sum())

        gnd[iv2plus] += 1
        self.assertEqual(2500,gnd.sum())
    
    def test_iadd_with_normalize_raises_warning(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        gnd = self.test_class(chroms)
        gnd.set_normalize(True)
        
        def my_func(ga):
            ga.set_normalize(True)
            ga[GenomicSegment("chrA",0,1000,"+")] += 5
      
        # manually reset registry before test
        plastid.util.services.exceptions.pl_once_registry = {}
        with warnings.catch_warnings(record=True) as warns:
            warnings.simplefilter("always")
            my_func(gnd)
        
        got_warning = False
        for w in warns:
            if "turning off normalization" in str(w.message):
                got_warning = True

        self.assertTrue(got_warning)
                
    def test_eq(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        gnd1 = self.test_class(chroms)
        gnd2 = self.test_class(chroms)

        self.assertEqual(gnd1,gnd2)
        
        # diff chrom nonzero positions
        gnd2[GenomicSegment("chrC",500,1000,"-")] = 200
        self.assertNotEqual(gnd1,gnd2)

        # same chroms, nonzero positions, and values
        gnd1[GenomicSegment("chrC",500,1000,"-")] = 200
        self.assertEqual(gnd1,gnd2)
        
        # same chroms and nonzero positions, diff values
        gnd1[GenomicSegment("chrC",500,1000,"-")] += 200
        self.assertNotEqual(gnd1,gnd2)
        
    def test_setters_and_getters(self):
        chroms = { "chrA" : 1000, "chrB" : 10000 }
        gnd = self.test_class(chroms)
        
        #generate random read data set
        for chr_name,chr_len in chroms.items():
            for strand in ("+","-"):
                num_iters = numpy.random.randint(50,100)#000)
                starts = numpy.random.randint(0,chr_len-50,size=num_iters)
                ends = [X + numpy.random.randint(0,50) for X in starts]
                vals = numpy.random.randint(0,100,size=num_iters)
                for i in range(num_iters):
                    iv = GenomicSegment(chr_name,int(starts[i]),int(ends[i]),strand)
                    gnd[iv] = vals[i]
        
        self.assertEquals(len(gnd),2*sum(chroms.values()),
                          "Chromosome lengths incorrect: %s vs %s" % (len(gnd),2*sum(chroms.values())))
        
        # auto-grow
        iv1 = GenomicSegment("chrA",10000,11000,"+")
        iv2 = GenomicSegment("chrA",10500,11000,"+")
        iv3 = GenomicSegment("chrA",int(5*1e5) + 10500,int(5*1e5) + 11000,"+")
        iv4 = GenomicSegment("chrB",int(5*1e5) + 10500,int(5*1e5) + 11000,"+")
        gnd[iv1] = 1
        self.assertEquals(sum(gnd[iv1]),1000,"Auto-grow failed during set")
        gnd[iv2] += 1
        self.assertEquals(sum(gnd[iv1]),1500,"+= failed during set")
        self.assertEquals(sum(gnd[iv2]),1000)

        self.assertGreater(gnd.lengths()["chrA"],1000,"Auto-grow failed chrA")

        gnd[iv3] += 1
        self.assertEqual(sum(gnd[iv3]),500,"+= failed during set")
        self.assertEqual(sum(gnd[iv4]),0,"Counts where not expected")

        # setters & getters
        iv1 = GenomicSegment("chrA",200,500,"+")
        oldvals = copy.deepcopy(gnd[iv1])
        gnd[iv1] = 5
        newvals = copy.deepcopy(gnd[iv1])
        self.assertTrue((oldvals != newvals).any(),"Set failed")
        self.assertEqual(newvals.sum(),5*len(newvals))
        self.assertEqual(newvals.sum(),5*len(iv1))

        newvals2 = copy.deepcopy(gnd[iv1])
        self.assertTrue((newvals2==newvals).all(),"Set failed")

        newvals = newvals2 = None
        gnd[iv1] = oldvals
        new_old = gnd[iv1]
        self.assertTrue((new_old==oldvals).all(),"Set valed")

        # scalar add
        gnd_plus5 = gnd + 5
        for iv in (iv1,iv2,iv3,iv4):
            self.assertTrue((gnd_plus5[iv] == (gnd[iv]+5)).all(),
                            "Scalar addition globally failed interval test")

        # FIXME- what is len of a sparse array?
        diff = abs(gnd_plus5.sum() - (gnd.sum() + 5*len(gnd)))
        self.assertLess(diff,self.tol,
                        "Error in genome-wide scalar addition (%s) exceeded tolerance %s. Raw vals: %s vs %s " % (diff,
                                                                                                                  self.tol,
                                                                                                                  gnd_plus5.sum(),
                                                                                                                  gnd.sum() + 5*len(gnd)))

        # scalar multiply
        gnd3 = gnd * 3
        for iv in (iv1,iv2,iv3,iv4):
            self.assertTrue((gnd3[iv] == (gnd[iv]*3)).all(),
                            "Scalar multip.lication globally failed interal test")

        quotient = 1.0 - (gnd3.sum() / 3*gnd.sum())
        self.assertLess(quotient,self.tol,
                        "Error in scalar multiplication (%s) exceeded tolerance %s" % (diff,self.tol))

        # genome-wide multiply
        gnd4 = gnd + gnd3
        gndmul = gnd3 * gnd4
        for iv in (iv1,iv2,iv3,iv4):
            self.assertTrue((gndmul[iv] - (gnd3[iv] * gnd4[iv]) <= self.tol).all(),
                            "Error in genome-wide multiply failed exceeded tolerance %s" % self.tol)

        # genome-wide add
        is_ok = True 
        for iv in (iv1,iv2,iv3):
            is_ok &= sum(gnd4[iv]) > 0
            self.assertTrue((gnd4[iv] == gnd[iv] + gnd3[iv]).all(),
                            "Error in genome-wide addition exceeded tolerance %s"  % self.tol)
            self.assertGreater(sum(gnd4[iv]),0)
    
    def test_nonzero(self):
        excluded_set = set([])
        for region in self.region_classes["repeat"] + self.region_classes["splice"]:
            excluded_set |= region.get_position_set()
        
        for k,v in self.gnds.items():
            nz_dict = v.nonzero()
            for strand in v.strands():
                strand_key = STRAND_KEYS[strand]
                expected = self.count_vecs["%s_%s" % (k,strand_key)].nonzero()[0]
                
                found = nz_dict["chrA"][strand]
                self.assertEqual(set(expected)-excluded_set,set(found)-excluded_set)

@attr(test="unit")
@attr(speed="slow")
class TestSparseGenomeArray(TestGenomeArray):
    """Test suite for |SparseGenomeArray|"""

    set_up   = False
    has_gnds = False

    def __init__(self,
                 methodName='runTest',
                 params=_SPARSE_GENOME_ARRAY_PARAMS,
                 test_folder=resource_filename("plastid","test/data/mini"),
                 tol=1e-8):
        """Initialize test case to run a single method. 
        We override this method to make sure expensive operations are only run when
        the first instance is made, and then stored in class attributes

        Parameters
        ----------
        methodName : str
            Name of method being run. Required by :py:class:`unittest.TestCase`
        
        params : dict
            Parameters specific to the set-up of test suites for specific
            AbstractgenomeArray subclasses. Don't change these
            
        test_folder : str or Resource
            Real or virtual location of folder of test data

        tol : float
            Tolerance for numerical differences between expected and observed
            values in the various tests
        """
        TestGenomeArray.__init__(self,methodName=methodName,params=params,test_folder=test_folder,tol=tol)



@attr(test="unit")
@attr(speed="slow")
class TestBigWigGenomeArray(AbstractGenomeArrayHelper):
    """Test suite for |SparseGenomeArray|"""

    set_up   = False
    has_gnds = False

    def __init__(self,
                 methodName='runTest',
                 params=_BIGWIG_GENOME_ARRAY_PARAMS,
                 test_folder=resource_filename("plastid","test/data/mini"),
                 tol=1e-3):
        """Initialize test case to run a single method. 
        We override this method to make sure expensive operations are only run when
        the first instance is made, and then stored in class attributes

        Parameters
        ----------
        methodName : str
            Name of method being run. Required by :py:class:`unittest.TestCase`
        
        params : dict
            Parameters specific to the set-up of test suites for specific
            AbstractgenomeArray subclasses. Don't change these
            
        test_folder : str or Resource
            Real or virtual location of folder of test data

        tol : float
            Tolerance for numerical differences between expected and observed
            values in the various tests
        """
        AbstractGenomeArrayHelper.__init__(self,methodName=methodName,params=params,test_folder=test_folder,tol=tol)
        if self.__class__.has_gnds == False:
            TestBigWigGenomeArray.gnds = TestBigWigGenomeArray.read_bigwig_files()
            self.__class__.has_gnds = True    
 
    @staticmethod
    def read_bigwig_files():
        """Read bigwig files into a dictionary
        """
        dtmp = {}
        for k in _SAMPLE_BASES:
            dtmp[k] = ga = BigWigGenomeArray(fill=0)
            
            fw = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_%s_fw.bw" % k)
            rc = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_%s_rc.bw" % k)

            ga.add_from_bigwig(fw, "+")
            ga.add_from_bigwig(rc, "-")
        
        return dtmp
            
#     def test_fill_value(self):
#         assert False
         
    def test_multiple_same_strand_sum(self):
        # should see sum double
        bigwigfile = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_center_12_fw.bw")
        bw = BigWigGenomeArray(fill=0)

        bw.add_from_bigwig(bigwigfile, "+")
        self.assertLessEqual(abs(bw.sum() - 4000),self.tol)

        bw.add_from_bigwig(bigwigfile, "+")
        self.assertLessEqual(abs(bw.sum() - 8000),self.tol)

        bw.add_from_bigwig(bigwigfile, "+")
        self.assertLessEqual(abs(bw.sum() - 12000),self.tol)

    def test_multiple_same_strand_fetch(self):
        bigwigfw = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_center_12_fw.bw")
        bigwigrc = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_center_12_rc.bw")
        wigfw = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_center_12_fw.wig")
        wigrc = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_center_12_rc.wig")
        
        bw = BigWigGenomeArray(fill=0)
        bw.add_from_bigwig(bigwigfw, "+")
        bw.add_from_bigwig(bigwigfw, "+")
        bw.add_from_bigwig(bigwigrc, "-")
        bw.add_from_bigwig(bigwigrc, "-")
        
        ga = GenomeArray(bw.lengths())
        ga.add_from_wiggle(open(wigfw),"+")
        ga.add_from_wiggle(open(wigrc),"-")
        
        for chrom, length in bw.lengths().items():
            for strand in bw.strands():
                seg = GenomicSegment(chrom,0,length,strand)
                maxdiff = abs(bw[seg] - 2*ga[seg]).max()
                msg = "Maximum difference for multiple_strand_fetch (%s) exceeds tolerance (%s)"% (maxdiff,self.tol)
                self.assertLessEqual(maxdiff, self.tol, msg)
     
    def test_to_genome_array(self):
        for test, orig in self.gnds.items():
            fw = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_%s_fw.wig" % test)
            rc = os.path.join(TestBigWigGenomeArray.test_folder,"wig","bw_%s_rc.wig" % test)            
            
            expected = GenomeArray()
            expected.add_from_wiggle(open(fw), "+")
            expected.add_from_wiggle(open(rc), "-")
            
            found = orig.to_genome_array()
            
            for chrom, length in expected.lengths().items():
                for strand in ("+","-"):
                    seg = GenomicSegment(chrom,0,length,strand)
                    diffvec = abs(orig[seg] - found[seg])
                    diffmax = diffvec.max() 
                    msg1 = "Maximum difference between exported GenomeArray and BigWigGenomeArray (%s) exceeds tolerance (%s) for test '%s' strand '%s'" % (diffmax,self.tol,test,strand)
                    self.assertLessEqual(diffmax,self.tol,msg1)
                    
                    
            for chrom, length in expected.lengths().items():
                for strand in ("+","-"):
                    seg = GenomicSegment(chrom,0,length,strand)
                    diffvec = abs(expected[seg] - found[seg])
                    diffmax = diffvec.max() 
                    msg1 = "Maximum difference between exported GenomeArray and wiggle-imported array (%s) exceeds tolerance (%s) for test '%s' strand '%s'" % (diffmax,self.tol,test,strand)
                    
                    self.assertLessEqual(diffmax,self.tol,msg1)
        
 
class FakeDict(object):
    """Creates a dictionary-like object that provies dictionary-like access
    to a BAMGenomeArray under various mapping rules, as if it were a collection
    of separate GenomeArrays. This is only a convenience class to allow us to
    re-use functions in the |AbstractGenomeArrayHelper| test suite in
    |TestBAMGenomeArray|
    """

    def __init__(self,bga,map_functions=_BAM_MAP_RULES):
        """Create a FakeDict

        Parameters
        ----------
        bga : |BAMGenomeArray|

        map_functions : dict
            Dictionary mapping descriptive names to mapping functions,
            such as those made by :py:func:`plastid.genomics.genome_array.FivePrimeMapFactory`
        """
        self.bga = bga
        self.map_functions = map_functions

    def __getitem__(self,key):
        self.bga.set_mapping(self.map_functions[key])
        return self.bga

    def items(self):
        for k in self.map_functions:
            yield (k,self[k])
    
    def values(self):
        # must use key, to trigger map setting in __getitem__
        for k in self.map_functions:
            yield self[k]
            
               
@attr(test="unit")
@attr(speed="slow")
class TestBAMGenomeArray(AbstractExportableGenomeArrayHelper):
    
    set_up   = False
    has_gnds = False

    def __init__(self,methodName='runTest',
                 params=_BAM_GENOME_ARRAY_PARAMS,
                 test_folder=resource_filename("plastid","test/data/mini"),
                 tol=1e-8):
        """Initialize test case to run a single method. 
        We override this method to make sure expensive operations are only run when
        the first instance is made, and then stored in class attributes

        Parameters
        ----------
        methodName : str
            Name of method being run. Required by py:class:`unittest.TestCase`

        params : dict
            Parameters specific to the set-up of test suites for specific
            AbstractgenomeArray subclasses. Don't change these
            
        test_folder : str or Resource
            Real or virtual location of folder of test data

        tol : float
            Tolerance for numerical differences between expected and observed
            values in the various tests
        """
        AbstractGenomeArrayHelper.__init__(self,methodName=methodName,params=params,test_folder=test_folder,tol=tol)
        if self.__class__.has_gnds == False:
            bga = BAMGenomeArray([pysam.Samfile(os.path.join(self.test_folder,_TEST_FILES["bam"]),"rb")])
            TestBAMGenomeArray.gnds = FakeDict(bga)
            TestBAMGenomeArray.bga = bga
            self.__class__.has_gnds = True

    def test_open_str_filename(self):
        z = BAMGenomeArray(os.path.join(self.test_folder,_TEST_FILES["bam"]))
        self.assertEqual(z.sum(),self.bga.sum())

    def test_open_multi_list(self):
        v = [os.path.join(self.test_folder,_TEST_FILES["bam"])]*2
        z = BAMGenomeArray(v)
        self.assertEqual(z.sum(),2*self.bga.sum())

    def test_open_multi_filename(self):
        f = os.path.join(self.test_folder,_TEST_FILES["bam"])
        z = BAMGenomeArray(f,f)
        self.assertEqual(z.sum(),2*self.bga.sum())

        
    def mutable_conversion_helper(self,new_class):
        """Helper function to test conversion of |BAMGenomeArray| to various |MutableAbstractGenomeArray| types

        Parameters
        ----------
        new_class : class
            Non-abstract subclass of |MutableAbstractGenomeArray|
        """
        ivplus = GenomicSegment("chrA",0,self.bga.lengths()["chrA"],"+")
        ivminus = GenomicSegment("chrA",0,self.bga.lengths()["chrA"],"-")
        for k,v in self.gnds.items():
            new_gnd = v.to_genome_array(array_type=new_class)

            for iv in (ivplus,ivminus):

                self.assertGreater(v[iv].sum(),0)
                self.assertGreater(new_gnd[iv].sum(),0)

                max_err = max(abs(v[iv]-new_gnd[iv]))
                err_message = "%s BAMGenomeArray conversion to %s error %s exceeded tolerance %s." % (k,new_class.__name__,max_err,self.tol)
                self.assertLess(max_err,self.tol,err_message)

    def test_to_genome_array(self):
        self.mutable_conversion_helper(GenomeArray)

    def test_to_sparse_genome_array(self):
        self.mutable_conversion_helper(SparseGenomeArray)

    def variablestep_and_bedgraph_export_helper(self,wiggle_type,export_function):
        # override function so we can test window size parameters in export
        for window_size in (1,2,5,10,25,100,500,1000,10000):
            AbstractGenomeArrayHelper.variablestep_and_bedgraph_export_helper(self,
                                                                              wiggle_type,
                                                                              export_function,
                                                                              input_class=GenomeArray,
                                                                              window_size=window_size)
    def test_add_remove_filter(self):
        # add a plus-strand filter and require minus strand regions be zero
        # then remove and watch it come back
        bga = self.bga
        
        def minus_only_filter(read):
            return read.is_reverse
        
        entire_iv_plus  = GenomicSegment("chrA",0,bga.lengths()["chrA"],"+")
        entire_iv_minus = GenomicSegment("chrA",0,bga.lengths()["chrA"],"-")
        
        # fetch counts & check
        pre_plus  = bga[entire_iv_plus]
        self.assertGreater(pre_plus.sum(),0)
        
        pre_minus = bga[entire_iv_minus]
        self.assertGreater(pre_minus.sum(),0)
        
        # add filter, re-fetch
        bga.add_filter("minus_only",minus_only_filter)
        post_plus  = bga[entire_iv_plus]
        post_minus = bga[entire_iv_minus]
        
        self.assertEqual(post_plus.sum(),0)
        self.assertFalse((post_plus==pre_plus).all())
        
        self.assertEqual(post_minus.sum(),pre_minus.sum())
        self.assertTrue((post_minus==pre_minus).all())
        
        # remove_filter, re_fetch
        bga.remove_filter("minus_only") 
        post_post_plus  = bga[entire_iv_plus]
        post_post_minus = bga[entire_iv_minus]

        self.assertEqual(post_post_plus.sum(),pre_plus.sum())
        self.assertTrue((post_post_plus==pre_plus).all())

        self.assertEqual(post_post_minus.sum(),pre_minus.sum())
        self.assertTrue((post_post_minus==pre_minus).all())

    

#===============================================================================
# INDEX: tools for generating test datasets with known results 
#===============================================================================

def _detect_or_create_folders(base_folder):
    """Creates and tests folder hierarchy needed for unit/integrative tests below.
    
    Parameters
    ----------
    base_folder : str
        path to base folder in which test data will be created
    """
    sub_folders = ["fasta","ebwt","count_vectors","align","wig","bed"]
    if not os.path.isdir(base_folder):
        os.mkdir(base_folder)
    for name in sub_folders:
        sf = os.path.join(base_folder,name)
        if not os.path.isdir(sf):
            os.mkdir(sf)
    
    # BED file for use later        
    fout = open(os.path.join(base_folder,_TEST_FILES["bed"]),"w") 
    fout.write(TEST_CHR_BED)
    fout.close()
    
    # .juncs file for tophat
    fout = open(os.path.join(base_folder,_TEST_FILES["juncs"]),"w")
    fout.write(TEST_CHR_JUNCS)
    fout.close()

def detect_base_folder(func):
    """Decorator function to ensure that folders required by functions below exist.
    For this decorator to work, the function it wraps MUST require base_folder
    as its first parameter.
    
    Parameters
    ----------
    func : Function
        Function to decorate
    
    Returns
    -------
    Function : wrapped function
    """
    @functools.wraps(func)
    def new_func(*args,**kwargs):
        _detect_or_create_folders(args[0])
        return func(*args,**kwargs)
    return new_func

@detect_base_folder
def create_synthetic_genome(base_folder):
    """Create a synthetic genome with unique and multimapping regions annotated
    as in TEST_CHR_BED
    
    Parameters
    ----------
    base_folder : str
        path to base folder in which test data will be created
        Genome will base saved as base_folder/fasta/chrA.fa
    
    Returns
    -------
    str : synthetic genome sequence
    """
    ustart  = 100
    ulength = 1000
    spacer_length = 100
    
    #uend = ustart + ulength
    #dstart = uend + spacer_length
    dlength = 500
    dspacer = 50
    
    #dstart = uend + spacer_length
    #dend   = dstart + dlength + dspacer + dlength
    
    #splice_start = dend + spacer_length
    splice_flank_length = 25
    splice_spacer = 75
    
    # generate sequence
    unique_region     = random_seq(ulength)
    duplicated_repeat = random_seq(dlength)
    duplicated_region = duplicated_repeat + ("N"*dspacer) + duplicated_repeat
    splice_region     = random_seq(splice_flank_length) + ("N"*splice_spacer) + random_seq(splice_flank_length)
    
    # write to FASTA
    total_genome = ("N"*ustart) + unique_region + "N"*spacer_length + duplicated_region + "N"*spacer_length + splice_region + "N"*100
    fh = open(os.path.join(base_folder,_TEST_FILES["genome"]),"w")
    fh.write(">%s\n%s\n" % ("chrA",total_genome))
    fh.close()
    return total_genome

@detect_base_folder
def create_bowtie_index(base_folder,bowtie_location=os.path.join(os.path.sep,"usr","local","bin")):
    """Build bowtie indices to enable alignments in bowtie and tophat
    against synthetic genome. Creates bowtie indices in base_folder/ebwt. 
    Requires a FASTA file of the synthetic genome in base_folder/fasta/chrA.fa,
    so run create_synthetic_genome() first

    Parameters
    ----------
    base_folder : str
        path to base folder in which test data will be created
    
    bowtie_location : str
        path to folder containing bowtie-build
    
    Returns
    -------
    int : exit status of bowtie-build
    """
    unspliced_args = [os.path.join(bowtie_location,"bowtie-build"),
                      os.path.join(base_folder,_TEST_FILES["genome"]),
                      os.path.join(base_folder,_TEST_FILES["bowtie_index"])]

    unspliced_exit = subprocess.call(unspliced_args)
    return unspliced_exit #| spliced_exit

def _ndarray_to_variable_step(count_vec,fh,name):
    """Write a numpy.ndarray to a variableStep wiggle file
    
    Parameters
    ----------
    count_vec : numpy.ndarray
        vector of counts
    
    fh : file-like
        open filehandle
    
    name : str
        Track name
    """
    fh.write("track type=wiggle_0 name=%s\n" % name)
    fh.write("variableStep chrom=chrA span=1\n")
    for i in count_vec.nonzero()[0]:
        val = count_vec[i]
        fh.write("%s\t%s\n" % (i+1,val))

def _ndarray_to_bedgraph(count_vec,fh,name):
    """Write a numpy.ndarray to a BEDGraph file
    
    Parameters
    ----------
    count_vec : numpy.ndarray
        vector of counts
    
    fh : file-like
        open filehandle
    
    name : str
        Track name
    """
    fh.write("track type=bedGraph name=%s\n" % name)
    last_val = count_vec[0]
    start_i  = 0
    for i,val in enumerate(count_vec):
        if val != last_val:
            fh.write("%s\t%s\t%s\t%s\n" % ("chrA",start_i,i,last_val))
            start_i  = i
            last_val = val
    fh.write("%s\t%s\t%s\t%s\n" % ("chrA",start_i,i+1,last_val))

@detect_base_folder
def generate_reads(base_folder,
                   reads_per_region=DEFAULT_READS_PER_REGION,
                   read_length=DEFAULT_READ_LENGTH):
    """Generates 30-nucleotide reads from a genome created by create_synthetic_genome,
    choosing from uniquely-mapping and multimapping regions annotated in 
    TEST_CHR_BED. 10000 reads are generated for each region type. Reads are
    returned in FASTA format. Also saves a numpy array of how many reads are expected
    to align to each nucleotide position in the synthetic genome, if reads are mapped
    at their 5' ends.

    Parameters
    ----------
    base_folder : str
        path to base folder in which test data will be created.
        Reads will base saved as base_folder/fasta/chrA_reads.fa.
        Count vectors will base saved as various text files in
        base_folder/count_vectors    

    reads_per_region : int
        Number of reads to generate in each region
    
    read_length : int
        Length of reads to generate
    
    Returns
    -------
    dict : dictionary of numpy.ndarrays corresponding to expected number of counts
           at each genomic position, under various read alignment mapping rules 
    """
    genome = SeqIO.to_dict(SeqIO.parse(open(os.path.join(base_folder,"fasta","chrA.fa")),
                                       "fasta"))
    
    len_A = len(genome["chrA"])
    count_vectors = { "fiveprime_0_fw"   : numpy.zeros(len_A).astype(int),
                      "fiveprime_15_fw"  : numpy.zeros(len_A).astype(int),
                      "threeprime_0_fw"  : numpy.zeros(len_A).astype(int),
                      "threeprime_15_fw" : numpy.zeros(len_A).astype(int),
                      "center_0_fw"      : numpy.zeros(len_A),
                      "center_12_fw"     : numpy.zeros(len_A),
                      
                      "fiveprime_0_rc"   : numpy.zeros(len_A).astype(int),
                      "fiveprime_15_rc"  : numpy.zeros(len_A).astype(int),
                      "threeprime_0_rc"  : numpy.zeros(len_A).astype(int),
                      "threeprime_15_rc" : numpy.zeros(len_A).astype(int),
                      "center_0_rc"      : numpy.zeros(len_A),
                      "center_12_rc"     : numpy.zeros(len_A),                     
                     }
    
    read_fh  = open(os.path.join(base_folder,_TEST_FILES["reads"]),"w")
    
    regions = filter(lambda x: "intron" not in x.get_name()\
                               and "entire" not in x.get_name(),
                               fetch_regions())
    
    for region in regions:
        strand_key = STRAND_KEYS[region.spanning_segment.strand]
        my_seq    = region.get_sequence(genome)
        
        # choose 5' read locations
        read_locs = numpy.random.randint(0,high=len(my_seq)-read_length+1,size=reads_per_region)
        
        # generate FASTA File
        # and save read positions to count vectors under various alignment mapping rules
        for n,loc in enumerate(read_locs):
            # write reads
            read_fh.write(">%s_%s\n%s\n" % (region.get_name(),n,my_seq[loc:loc+read_length]))
            _,position,_ = region.get_genomic_coordinate(loc)

            # populate 5' and 3' mapped count vectors
            for offset in (0,15):
                _,position,_ = region.get_genomic_coordinate(loc + offset)
                count_vectors["fiveprime_%s_%s" % (offset,strand_key)][position] += 1
                _,position,_ = region.get_genomic_coordinate(loc + read_length - offset - 1 )
                count_vectors["threeprime_%s_%s" % (offset,strand_key)][position] += 1
            
            # populate center-mapped count vectors
            read_positions = region.get_subchain(loc,loc+read_length).get_position_list()
            assert len(read_positions) == read_length
            for pos in read_positions:
                count_vectors["center_0_%s" % strand_key][pos] += 1.0/len(read_positions)
            assert len(read_positions[12:-12]) == read_length - 24
            for pos in read_positions[12:-12]:
                count_vectors["center_12_%s" % strand_key][pos] += 1.0/(len(read_positions)-24)

    for k,v in count_vectors.items():
        numpy.savetxt(os.path.join(base_folder,"count_vectors","%s.txt" % k),v)

    read_fh.close()
    
    # export 5' mapped BEDGraph files
    bedgraph_fw = open(os.path.join(base_folder,_TEST_FILES["bedgraph_fw"]),"w")
    bedgraph_rc = open(os.path.join(base_folder,_TEST_FILES["bedgraph_rc"]),"w")
    
    _ndarray_to_bedgraph(count_vectors["fiveprime_0_fw"],bedgraph_fw,base_folder)
    _ndarray_to_bedgraph(count_vectors["fiveprime_0_rc"],bedgraph_rc,base_folder)
    
    bedgraph_fw.close()
    bedgraph_rc.close()

    # export 5' mapped variableStep wiggle files
    vs_fw = open(os.path.join(base_folder,_TEST_FILES["variable_step_fw"]),"w")
    vs_rc = open(os.path.join(base_folder,_TEST_FILES["variable_step_rc"]),"w")

    _ndarray_to_variable_step(count_vectors["fiveprime_0_fw"],vs_fw,base_folder)
    _ndarray_to_variable_step(count_vectors["fiveprime_0_rc"],vs_rc,base_folder)

    vs_fw.close()
    vs_rc.close()
    
    return count_vectors

@detect_base_folder
def perform_alignments(base_folder,
                       bowtie_location=os.path.join(os.path.sep,"usr","local","bin"),
                       tophat_location=os.path.join(os.path.sep,"usr","local","bin"),
                       samtools_location=os.path.join(os.path.sep,"usr","local","bin"),
                       ):
    """Perform alignments of generated reads against synthetic genome,
    in both Tophat and Bowtie so that both BAM and bowtie input may be tested.
    Note: spliced reads will not align in bowtie.
    
    Read alignments will be placed in base_folder/align
    
    Requires a bowtie index of the synthetic genome in base_folder/ebwt, and
    syntheti reads to align in base_folder/fasta/chrA_reads.fa 
    so run build_bowtie_index() and generate_reads() first.

    Parameters
    ----------
    base_folder : str
        path to base folder in which test data will be created.
    
    bowtie_location : str
        path to folder containing bowtie executable

    tophat_location : str
        path to folder containing tophat executable

    samtools_location : str
        path to folder containing samtools executable
    
    Returns
    -------
    int : ORed exit status of bowtie and tophat (0 if both successful; 1 otherwise)
    """
    # align with no mismatches, choosing 1 random alignment from repeat regions
    bowtie_args = [os.path.join(bowtie_location,"bowtie"),
                   "-v0","-k1","--best","-f",
                   "--un",os.path.join(base_folder,"align","chrA_unspliced_unaligned.fa"),
                   os.path.join(base_folder,_TEST_FILES["bowtie_index"]),
                   os.path.join(base_folder,_TEST_FILES["reads"]),
                   os.path.join(base_folder,_TEST_FILES["bowtie"])
                   ]
    
    # align with no mismatches, choosing 1 random alignment from repeat regions
    tophat_args = [os.path.join(bowtie_location,"tophat"),
                   "--bowtie1",
                   "--read-mismatches=0",
                   "--min-intron-length=20",
                   "--library-type=fr-firststrand",
                   "--raw-juncs",
                   os.path.join(base_folder,_TEST_FILES["juncs"]),
                   "--no-novel-juncs",
                   "-o",
                   os.path.join(base_folder,"align","tophat"),
                   os.path.join(base_folder,_TEST_FILES["bowtie_index"]),
                   os.path.join(base_folder,_TEST_FILES["reads"]),
                   ]

    samtools_multi_args = [os.path.join(samtools_location,"samtools"),
                           "view","-F", "256",
                           os.path.join(base_folder,"align","tophat","accepted_hits.bam"),
                           "-b",
                           "-o",
                           os.path.join(base_folder,"align","chrA_tophat.bam")
                          ]

    samtools_index_args = [os.path.join(samtools_location,"samtools"),
                           "index",
                           os.path.join(base_folder,"align","chrA_tophat.bam")
                          ]
    
    bowtie_exit   = subprocess.call(bowtie_args)
    tophat_exit   = subprocess.call(tophat_args)
    
    samtools_exit_1 = subprocess.call(samtools_multi_args)
    samtools_exit_2 = subprocess.call(samtools_index_args)
    
    return bowtie_exit | tophat_exit | samtools_exit_1 | samtools_exit_2

def create_dataset(base_folder,
                   bowtie_location=os.path.join(os.path.sep,"usr","local","bin"),
                   tophat_location=os.path.join(os.path.sep,"usr","local","bin"),
                   samtools_location=os.path.join(os.path.sep,"usr","local","bin"),):
    """Create a ground-truth dataset for testing |GenomeArray|
    
    This dataset includes a synthetic genome of random sequence, containing
    uniquely-mapping and multimapping regions; short sequence reads generated
    from these regions; alignments of these reads made in bowtie and tophat;
    and saved numpy tables indicating how many counts should appear at each
    genomic position under various mapping rules and offsets.
    
    Parameters
    ----------
    base_folder : str
        path to base folder in which test data will be created.
    
    bowtie_location : str
        path to folder containing bowtie executable

    tophat_location : str
        path to folder containing tophat executable

    samtools_location : str
        path to folder containing samtools executable
    
    Returns
    -------
    dict: dictionary containing the genome sequence, count vectors,
          and alignment status
    """
    dtmp = {}
    dtmp["base_folder"] = base_folder
    dtmp["genome_str"]  = create_synthetic_genome(base_folder)
    dtmp["aligned"]     = False
    if create_bowtie_index(base_folder,bowtie_location=bowtie_location) == 0:
        dtmp.update(generate_reads(base_folder))
    
        if perform_alignments(base_folder,
                              bowtie_location=bowtie_location,
                              tophat_location=tophat_location,
                              samtools_location=samtools_location) == 0:
            dtmp["aligned"] = True
            return dtmp
    
    return dtmp
