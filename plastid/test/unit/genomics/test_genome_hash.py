#!/usr/bin/env python
"""Tests for data structures defined in :py:mod:`plastid.genomics.genome_hash`
"""
import unittest
from random import shuffle
from pkg_resources import resource_filename, cleanup_resources
from nose.plugins.attrib import attr

from plastid.genomics.roitools import Transcript
from plastid.readers.bed import BED_Reader
from plastid.genomics.genome_hash import GenomeHash, BigBedGenomeHash, TabixGenomeHash
from plastid.util.services.decorators import skip_if_abstract
from plastid.util.io.filters import CommentReader
from plastid.readers.bigbed import BigBedReader

from plastid.test.ref_files import REF_FILES
import warnings
warnings.simplefilter("ignore",DeprecationWarning)

#===============================================================================
# INDEX: Test case for GenomeHash
#===============================================================================    

def tearDownModule():
    """Remove test dataset files after unit tests are complete"""
    cleanup_resources()

def _name_sort(x):
    return str(x)

class AbstractGenomeHashHelper(unittest.TestCase):
   
    @skip_if_abstract
    def test_getitem_genomicsegment(self):
#        """test getitem
#        
#        1.  Make sure each feature can fetch its own subregion from its own neighborhood
#        
#        2.  Make sure each feature cannot fetch its own antisense subregion
#        
#        3.  Make sure that a 5' UTR cannot fetch its corresponding CDS
#        """             
        # make sure we can fetch each transcript's own CDS
        u = 0
        for txid,tx in self.tx_dict.items():
            ol_features = list(self.cds_hash[tx.spanning_segment])
            self.assertIn(txid,(X.get_name() for X in ol_features),
                          msg="%s failed to fetch its own CDS on same strand. Got %s" % (txid,ol_features))
        
            # make sure 5' UTR does not fetch CDS
            utr5 = tx.get_utr5()
            if len(utr5) > 0:
                u += 1
                ol_features = list(self.cds_hash[utr5.spanning_segment])
                self.assertNotIn(txid,(X.get_name() for X in ol_features),
                                 msg="%s 5' UTR fetched CDS on same strand" % txid)
    
                ol_features = list(self.as_cds_hash[utr5.spanning_segment])
                self.assertNotIn(txid,(X.get_name() for X in ol_features),
                                 msg="%s 5' UTR fetched CDS on opposite strand" % txid)
            
            self.assertGreater(u,0)

        # make sure we don't fetch each transcript's own antisense CDS
        # on opposite strand
        for txid,tx in self.tx_dict.items():
            ol_features = list(self.as_cds_hash[tx.spanning_segment])
            self.assertNotIn(txid,(X.get_name() for X in ol_features),
                          msg="%s fetched its own name on wrong strand!" % txid)        



    @skip_if_abstract
    def test_getitem_segmentchain(self):
#        """test getitem
#        
#        1.  Make sure each feature can fetch its own subregion from its own neighborhood
#        
#        2.  Make sure each feature cannot fetch its own antisense subregion
#        
#        3.  Make sure that a 5' UTR cannot fetch its corresponding CDS
#        """             
        # make sure we can fetch each transcript's own CDS
        u = 0
        for txid,tx in self.tx_dict.items():
            ol_features = list(self.cds_hash[tx])
            self.assertIn(txid,(X.get_name() for X in ol_features),
                          msg="%s failed to fetch its own CDS on same strand. Got %s" % (txid,ol_features))
        
            # make sure 5' UTR does not fetch CDS
            utr5 = tx.get_utr5()
            if len(utr5) > 0:
                u += 1
                ol_features = list(self.cds_hash[utr5])
                self.assertNotIn(txid,(X.get_name() for X in ol_features),
                                 msg="%s 5' UTR fetched CDS on same strand" % txid)
    
                ol_features = list(self.as_cds_hash[utr5])
                self.assertNotIn(txid,(X.get_name() for X in ol_features),
                                 msg="%s 5' UTR fetched CDS on opposite strand" % txid)
            
            self.assertGreater(u,0)

        # make sure we don't fetch each transcript's own antisense CDS
        # on opposite strand
        for txid,tx in self.tx_dict.items():
            ol_features = list(self.as_cds_hash[tx])
            self.assertNotIn(txid,(X.get_name() for X in ol_features),
                          msg="%s fetched its own name on wrong strand!" % txid)        

   
    @skip_if_abstract
    def test_overlap_stranded(self):
#        """test get_overlapping_features with strand
#        
#        1.  Make sure each feature can fetch its own subregion from its own neighborhood
#        
#        2.  Make sure each feature cannot fetch its own antisense subregion
#        
#        3.  Make sure that a 5' UTR cannot fetch its corresponding CDS
#        """             
        # make sure we can fetch each transcript's own CDS
        u = 0
        for txid,tx in self.tx_dict.items():
            ol_features = list(self.cds_hash.get_overlapping_features(tx,stranded=True))
            self.assertIn(txid,(X.get_name() for X in ol_features),
                          msg="%s failed to fetch its own CDS on same strand. Got %s" % (txid,ol_features))
        
            # make sure 5' UTR does not fetch CDS
            utr5 = tx.get_utr5()
            if len(utr5) > 0:
                u += 1
                ol_features = list(self.cds_hash.get_overlapping_features(utr5,stranded=False))
                self.assertNotIn(txid,(X.get_name() for X in ol_features),
                                 msg="%s 5' UTR fetched CDS on same strand" % txid)
    
                ol_features = list(self.as_cds_hash.get_overlapping_features(utr5,stranded=False))
                self.assertNotIn(txid,(X.get_name() for X in ol_features),
                                 msg="%s 5' UTR fetched CDS on opposite strand" % txid)
            
            self.assertGreater(u,0)

        # make sure we don't fetch each transcript's own antisense CDS
        # on opposite strand
        for txid,tx in self.tx_dict.items():
            ol_features = list(self.as_cds_hash.get_overlapping_features(tx,stranded=True))
            self.assertNotIn(txid,(X.get_name() for X in ol_features),
                          msg="%s fetched its own name on wrong strand!" % txid)    
    
    @skip_if_abstract
    def test_overlap_unstranded(self):
#        """test get_overlapping_features unstranded
#        
#        1.  Make sure each feature can fetch its own subregion from its own neighborhood
#        
#        2.  Make sure each feature can fetch its own antisense subregion
#        
#        3.  Make sure that a 5' UTR cannot fetch its corresponding CDS
#        """        
        u = 0
        for txid,tx in self.tx_dict.items():
            # make sure we can fetch each transcript's own CDS from both sense and antisense
            ol_features = list(self.cds_hash.get_overlapping_features(tx,stranded=False))
            self.assertIn(txid,(X.get_name() for X in ol_features),
                          msg="%s failed to fetch its own name on same strand." % txid)

            ol_features = list(self.as_cds_hash.get_overlapping_features(tx,stranded=False))
            self.assertIn(txid,(X.get_name() for X in ol_features),
                          msg="%s failed to fetch its own name on opposite strand." % txid)

            utr5 = tx.get_utr5()
            if len(utr5) > 0:
                u += 1
                ol_features = list(self.cds_hash.get_overlapping_features(utr5,stranded=False))
                self.assertNotIn(txid,(X.get_name() for X in ol_features),
                                 msg="%s 5' UTR fetched CDS on same strand" % txid)
    
                ol_features = list(self.as_cds_hash.get_overlapping_features(utr5,stranded=False))
                self.assertNotIn(txid,(X.get_name() for X in ol_features),
                                 msg="%s 5' UTR fetched CDS on opposite strand" % txid)
            
            self.assertGreater(u,0)

@attr(test="unit")
class TestGenomeHash(AbstractGenomeHashHelper):
    """Test case for :py:class:`GenomeHash`"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data for `TestGenomeHash`"""
        cls.binsize = 10000

        cls.transcripts    = list(BED_Reader(CommentReader(open(REF_FILES["100transcripts_bed"])),return_type=Transcript))
        cls.coding_regions = list(BED_Reader(CommentReader(open(REF_FILES["100cds_bed"])),return_type=Transcript))
        cls.coding_antisense = list(BED_Reader(CommentReader(open(REF_FILES["100cds_antisense_bed"])),return_type=Transcript))

        cls.tx_dict  = { X.get_name() : X for X in cls.transcripts }
        cls.cds_dict = { X.get_name() : X for X in cls.coding_regions }
        cls.as_cds_dict = { X.get_name() : X for X in cls.coding_antisense }

        cls.tx_hash     = GenomeHash(cls.tx_dict,do_copy=False,binsize=cls.binsize)
        cls.cds_hash    = GenomeHash(cls.cds_dict,do_copy=False,binsize=cls.binsize)
        cls.as_cds_hash = GenomeHash(cls.as_cds_dict,do_copy=False,binsize=cls.binsize)
        
        cls.shuffled_indices = list(range(len(cls.transcripts)))
        shuffle(cls.shuffled_indices)
    
    def test_genomehash_create_from_dict(self):
#        """Test creation `GenomeHash`es from dictionaries without loss of information"""
        gh = GenomeHash(self.tx_dict,do_copy=True)
        
        found    = sorted(list(gh.feature_dict.values()),key=_name_sort)
        expected = sorted(list(self.tx_dict.values()),key=_name_sort)
        self.assertEquals(len(expected),len(found),"Features lost in creation of GenomeHash from dict. Expected %s, found %s." % (len(expected),len(found)))
        self.assertEquals(sorted(list(gh.feature_dict.values()),key=_name_sort),
                          sorted(list(self.tx_dict.values()),key=_name_sort),
                          "Features lost in creation of GenomeHash from dict")
        
    def test_genomehash_create_from_list(self):
#        """Test creation `GenomeHash`es from lists without loss of information"""
        gh = GenomeHash(self.transcripts,do_copy=True)
        found    = sorted(list(gh.feature_dict.values()),key=_name_sort)
        expected = sorted(list(self.tx_dict.values()),key=_name_sort)
        self.assertEquals(len(expected),len(found),"Features lost in creation of GenomeHash from list. Expected %s, found %s." % (len(expected),len(found)))
        self.assertEquals(expected,
                          found,
                          "Features lost in creation of GenomeHash from list:\nFirst:\n%s\nSecond:\n%s" % (expected,found))
    
    def test_get_hash_bins(self):
        for ivc in self.transcripts:
            bin_start = ivc.spanning_segment.start // self.binsize
            bin_end   = ivc.spanning_segment.end // self.binsize
            my_range = set(range(bin_start,bin_end+1))
            self.assertTrue(len(my_range)>0)
            self.assertEquals(my_range,set(self.tx_hash._get_hash_bins(ivc)))
    
    def test_genomehash_update_from_dict(self):
#        """Test addition to a `GenomeHash` from a dictionary without loss of features
#        
#        1. Add features to an empty dictionary
#        
#        2. Add features to a non-empty dictionary
#        """
        tuple_sort = lambda x: _name_sort(x[1])
        gh = GenomeHash({},do_copy=True)
        tx_items  = sorted(list(self.tx_dict.items()),key=tuple_sort)
        tx_values = [X[1] for X in tx_items]
        dict1 = dict(tx_items[:50])
        dict2 = dict(tx_items[50:])
        self.assertEqual(len(dict1),50)
        self.assertEqual(len(dict2),50)

        # check values, as opposed to items, because keys of feature_dict
        # are unique numerical IDs as opposed to ames
        gh.update(dict1)
        self.assertEquals(sorted(list(gh.feature_dict.values()),key=_name_sort),
                          tx_values[:50],
                          "Features lost in first update of empty GenomeHash from dict")
    
        gh.update(dict2)
        self.assertEquals(sorted(list(gh.feature_dict.values()),key=_name_sort),
                          tx_values,
                          "Features lost in second update of non-empty GenomeHash from dict")
    
    def test_genomehash_update_from_list(self):
#        """Test addition to a `GenomeHash` from a list without loss of features
#        
#        1. Add features to an empty dictionary
#        
#        2. Add features to a non-empty dictionary
#        """
        gh = GenomeHash({},do_copy=True)
        tx_list = sorted(list(self.tx_dict.values()),key=_name_sort)
        self.assertGreater(len(tx_list),0)

        gh.update(tx_list[:50])
        self.assertEquals(sorted(list(gh.feature_dict.values()),key=_name_sort),
                          tx_list[:50],
                          "Features lost in update of empty GenomeHash from list")
    
        gh.update(tx_list[50:])
        self.assertEquals(sorted(list(gh.feature_dict.values()),key=_name_sort),
                          tx_list,
                          "Features lost in update of non-empty GenomeHash from list")

    def test_nearby_feature_names_stranded(self):
#        """Test fetching of nearby feature names
#        
#        1.  Make sure each feature can fetch its own name from its own neighborhood
#        
#        2.  Make sure each feature cannot fetch its own name from antisense neighorhood
#        
#        2.  Make sure transcripts that are not on the same chromosome, not on
#            the same strand, or that are too far apart on the same chromosome
#            and strand, do not fetch each other.
#        """
        # make sure we can fetch each transcript's own CDS
        for txid,tx in self.tx_dict.items():
            nearby_names = self.cds_hash.get_nearby_feature_names(tx,stranded=True)
            self.assertIn(txid,nearby_names,
                          msg="%s failed to fetch its own name" % txid)

        # make sure we don't fetch each transcript's own antisense CDS
        # on opposite strand
        for txid,tx in self.tx_dict.items():
            nearby_names = self.as_cds_hash.get_nearby_feature_names(tx,
                                                                     stranded=True)
            self.assertNotIn(txid,nearby_names,
                          msg="%s fetched its own name on wrong strand!" % txid)            
            
        # make sure that at least one transcript doesn't fetch
        # another's CDS name if we know they are not:
        #     -on same chromosome
        #     -on same strand
        #     -too far on same chromosme and strand
        c = 0 # count how many times we would expect to fetch self
        for idx, (txid,tx) in zip(self.shuffled_indices,sorted(self.tx_dict.items())):
            fake_tx = self.tx_dict[sorted(self.tx_dict.keys())[idx]]
            if fake_tx.spanning_segment.chrom != tx.spanning_segment.chrom or\
               fake_tx.spanning_segment.strand != tx.spanning_segment.strand or\
                (abs(fake_tx.spanning_segment.start - tx.spanning_segment.start) > self.binsize and\
                 abs(fake_tx.spanning_segment.end - tx.spanning_segment.end) > self.binsize):
                nearby_names = self.cds_hash.get_nearby_feature_names(fake_tx)
                self.assertNotIn(txid,nearby_names)
            else:
                c += 1
        
        # make sure negative test was not trivial
        # i.e. that we tested at least some set of mismatched transcripts
        self.assertNotEqual(c,len(self.tx_dict))

    def test_nearby_feature_names_unstranded(self):
#        """Test fetching of nearby feature names, disregarding strand
#        
#        1.  Make sure each feature can fetch its own name from its own sense neighborhood
#        
#        2.  Make sure each feature can fetch its own name from antisense neighorhood
#        
#        2.  Make sure transcripts that are not on the same chromosome, or that are too
#            far apart on the same chromosome, do not fetch each other.
#        """ 
        # make sure we can fetch each transcript's own CDS from both sense asd antisense
        for txid,tx in self.tx_dict.items():
            nearby_names = self.cds_hash.get_nearby_feature_names(tx,stranded=False)
            self.assertIn(txid,nearby_names,
                          msg="%s failed to fetch its own name on same strand: %s" % (txid,
                                                                     ",".join(nearby_names)))

            nearby_names = self.as_cds_hash.get_nearby_feature_names(tx,stranded=False)
            self.assertIn(txid,nearby_names,
                          msg="%s failed to fetch its own name on opposite strand: %s" % (txid,
                                                                     ",".join(nearby_names)))

        # make sure that at least one transcript doesn't fetch
        # another's CDS name if we know they are not:
        #     -on same chromosome
        #     -too far on same chromosme and strand
        c = 0 # count how many times we would expect to fetch self
        for idx, (txid,tx) in zip(self.shuffled_indices,sorted(self.tx_dict.items())):
            fake_tx = self.tx_dict[sorted(self.tx_dict.keys())[idx]]
            if fake_tx.spanning_segment.chrom != tx.spanning_segment.chrom or\
                (abs(fake_tx.spanning_segment.start - tx.spanning_segment.start) > self.binsize and\
                 abs(fake_tx.spanning_segment.end - tx.spanning_segment.end) > self.binsize):
                nearby_names = self.cds_hash.get_nearby_feature_names(fake_tx)
                self.assertNotIn(txid,nearby_names)
            else:
                c += 1
        
        # make sure negative test was not trivial
        # i.e. that we tested at least some set of mismatched transcripts
        self.assertNotEqual(c,len(self.tx_dict))            
            
    def test_nearby_features_stranded(self):
#        """Test fetching of nearby features
#        
#        1.  Make sure each feature can fetch its own subregion from its own neighborhood
#        
#        2.  Make sure each feature cannot fetch its own antisense subregion
#        
#        3.  Make sure fatures that are not on the same chromosome, not on
#            the same strand, or that are too far apart on the same chromosome
#            and strand, do not fetch each other.
#        """
        # make sure we can fetch each transcript's own CDS
        for txid,tx in self.tx_dict.items():
            nearby_features = self.cds_hash.get_nearby_features(tx)
            self.assertIn(txid,(X.get_name() for X in nearby_features),
                          msg="%s failed to fetch its own CDS on correct strand" % txid)


        # make sure we don't fetch each transcript's own antisense CDS
        # on opposite strand
        for txid,tx in self.tx_dict.items():
            nearby_features = self.as_cds_hash.get_nearby_features(tx,stranded=True)
            self.assertNotIn(txid,(X.get_name() for X in nearby_features),
                          msg="%s fetched its own name on wrong strand!" % txid)    
            
        # make sure that at least one transcript doesn't fetch
        # another's CDS if we know they are not:
        #     -on same chromosome
        #     -on same strand
        #     -too far on same chromosme and strand
        c = 0 # count how many times we would expect to fetch self
        for idx, (txid,tx) in zip(self.shuffled_indices,sorted(self.tx_dict.items())):
            fake_tx = self.tx_dict[sorted(self.tx_dict.keys())[idx]]
            if fake_tx.spanning_segment.chrom != tx.spanning_segment.chrom or\
               fake_tx.spanning_segment.strand != tx.spanning_segment.strand or\
                (abs(fake_tx.spanning_segment.start - tx.spanning_segment.start) > self.binsize and\
                 abs(fake_tx.spanning_segment.end - tx.spanning_segment.end) > self.binsize):
                nearby_features = self.cds_hash.get_nearby_features(fake_tx)
                self.assertNotIn(txid,(X.get_name() for X in nearby_features))
            else:
                c += 1

        # make sure negative test was not trivial
        # i.e. that we tested at least some set of mismatched transcripts
        self.assertNotEqual(c,len(self.tx_dict))

    def test_nearby_feature_unstranded(self):
#        """Test fetching of nearby features, disregarding strand
#        
#        1.  Make sure each feature can fetch its own subregion from its own neighborhood
#        
#        2.  Make sure each feature can fetch its own antisense subregion
#        
#        3.  Make sure fatures that are not on the same chromosome, or that are too far
#            apart on the same chromosome, do not fetch each other.
#        """        
        # make sure we can fetch each transcript's own CDS from both sense asd antisense
        for txid,tx in self.tx_dict.items():
            nearby_features = self.cds_hash.get_nearby_features(tx,stranded=False)
            self.assertIn(txid,(X.get_name() for X in nearby_features),
                          msg="%s failed to fetch its own name on same strand." % txid)

            nearby_features = self.as_cds_hash.get_nearby_features(tx,stranded=False)
            self.assertIn(txid,(X.get_name() for X in nearby_features),
                          msg="%s failed to fetch its own name on opposite strand." % txid)


        # make sure that at least one transcript doesn't fetch
        # another's CDS if we know they are not:
        #     -on same chromosome
        #     -too far on same chromosme and strand
        c = 0 # count how many times we would expect to fetch self
        for idx, (txid,tx) in zip(self.shuffled_indices,sorted(self.tx_dict.items())):
            fake_tx = self.tx_dict[sorted(self.tx_dict.keys())[idx]]
            if fake_tx.spanning_segment.chrom != tx.spanning_segment.chrom or\
                (abs(fake_tx.spanning_segment.start - tx.spanning_segment.start) > self.binsize and\
                 abs(fake_tx.spanning_segment.end - tx.spanning_segment.end) > self.binsize):
                nearby_features = self.cds_hash.get_nearby_features(fake_tx,stranded=False)
                self.assertNotIn(txid,(X.get_name() for X in nearby_features))
            else:
                c += 1

@attr(test="unit")    
class TestBigBedGenomeHash(AbstractGenomeHashHelper):

    @classmethod
    def setUpClass(cls):
        """Set up test data for `TestGenomeHash`"""
        cls.binsize = 10000

        cls.tx_bbfile     = REF_FILES["100transcripts_bigbed"]
        cls.cds_bbfile    = REF_FILES["100cds_bigbed"]
        cls.as_cds_bbfile = REF_FILES["100cds_antisense_bigbed"]
        
        cls.tx_hash     = BigBedGenomeHash(cls.tx_bbfile)
        cls.cds_hash    = BigBedGenomeHash(cls.cds_bbfile)
        cls.as_cds_hash = BigBedGenomeHash(cls.as_cds_bbfile)
        
        cls.transcripts      = list(BigBedReader(cls.tx_bbfile,return_type=Transcript))
        cls.coding_regions   = list(BigBedReader(cls.cds_bbfile))
        cls.shuffled_indices = list(range(len(cls.transcripts)))
        
        cls.tx_dict  = { X.get_name() : X for X in cls.transcripts }
        cls.cds_dict = { X.get_name() : X for X in cls.coding_regions }
        shuffle(cls.shuffled_indices)

    def test_works_with_str_filename(self):
        one_hash = BigBedGenomeHash(self.tx_bbfile)
        for tx in self.transcripts:
            expected = self.tx_hash[tx]
            found = one_hash[tx]
            self.assertEqual(expected, found)
                    
    def test_works_with_multiple_files(self):
        # check to make sure all files are queried by comparing results
        # from a single hash to results from a double hash
        double_hash = BigBedGenomeHash(self.tx_bbfile,self.cds_bbfile)
        for tx in self.transcripts:
            expected = self.tx_hash[tx] + self.cds_hash[tx]
            found = double_hash[tx]
            self.assertEqual(expected, found)


@attr(test="unit")    
class TestTabixGenomeHash(AbstractGenomeHashHelper):
 
    @classmethod
    def setUpClass(cls):
        """Set up test data for `TestTabixGenomeHash`"""
 
        cls.tx_file     = REF_FILES["100transcripts_bed_tabix"]
        cls.cds_file    = REF_FILES["100cds_bed_tabix"]
        cls.as_cds_file = REF_FILES["100as_cds_bed_tabix"]

        cls.tx_hash     = TabixGenomeHash(cls.tx_file,data_format="BED")
        cls.cds_hash    = TabixGenomeHash(cls.cds_file,data_format="BED")
        cls.as_cds_hash = TabixGenomeHash(cls.as_cds_file,data_format="BED")

        # use BigBeds as reference        
        # TODO: change to Transcript objects
        cls.tx_bbfile     = REF_FILES["100transcripts_bigbed"]
        cls.cds_bbfile    = REF_FILES["100cds_bigbed"]
        cls.as_cds_bbfile = REF_FILES["100cds_antisense_bigbed"]
         
        cls.transcripts      = list(BigBedReader(cls.tx_bbfile,return_type=Transcript))
        cls.coding_regions   = list(BigBedReader(cls.cds_bbfile))
        cls.shuffled_indices = list(range(len(cls.transcripts)))
         
        cls.tx_dict  = { X.get_name() : X for X in cls.transcripts }
        cls.cds_dict = { X.get_name() : X for X in cls.coding_regions }
        shuffle(cls.shuffled_indices)
    
    def test_unknown_reader_raises_error(self):
        self.assertRaises(KeyError,TabixGenomeHash, self.tx_file, data_format="GARBAGE")
    
    def test_works_with_str_filename(self):
        one_hash = TabixGenomeHash(self.tx_file,data_format="BED")
        for tx in self.tx_dict.values():
            found = one_hash[tx]
            expected = self.tx_hash[tx]
            self.assertEqual(expected,found)
                    
    def test_works_with_multiple_files(self):
        # check to make sure all files are queried by comparing results
        # from a single hash to results from a double hash
        double_hash = TabixGenomeHash(self.tx_file,self.cds_file,data_format="BED")
        for tx in self.tx_dict.values():
            df = double_hash[tx]
            expected = self.tx_hash[tx] + self.cds_hash[tx]
            self.assertEqual(len(df),len(expected))

    def test_works_with_multiple_files_list(self):
        # check to make sure all files are queried by comparing results
        # from a single hash to results from a double hash
        double_hash = TabixGenomeHash([self.tx_file,self.cds_file],data_format="BED")
        for tx in self.tx_dict.values():
            df = double_hash[tx]
            expected = self.tx_hash[tx] + self.cds_hash[tx]
            self.assertEqual(len(df),len(expected))
