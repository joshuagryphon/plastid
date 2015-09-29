#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.readers.bigbed`

Notes
-----
Several of these tests are tested against |GenomeHash|, and so will fail if 
|GenomeHash| is malfunctioning
"""

import unittest
import copy
import warnings
from random import shuffle
from pkg_resources import resource_filename, cleanup_resources
from nose.plugins.attrib import attr
from collections import OrderedDict
from plastid.genomics.roitools import SegmentChain, GenomicSegment, Transcript
from plastid.genomics.genome_hash import GenomeHash
from plastid.readers.bed import BED_to_Transcripts, BED_to_SegmentChain, BED_Reader
from plastid.readers.bigbed import BigBedReader, RTree, RTreeLeafFactory

warnings.simplefilter("ignore",DeprecationWarning)

#===============================================================================
# INDEX: helper functions
#===============================================================================

def tearDownModule():
    """Remove test dataset files after unit tests are complete"""
    cleanup_resources()

def transcript_identical(ivc1,ivc2):
    """Test for identity between positions of two Transcripts"""
    position_test = ivc1.get_position_set() == ivc2.get_position_set()
    strand_test   = ivc1.spanning_segment.strand == ivc2.spanning_segment.strand
    chrom_test    = ivc1.spanning_segment.chrom == ivc2.spanning_segment.chrom

    start_test = (ivc1.cds_start is None and ivc2.cds_start is None) or\
                 (ivc1.cds_start == ivc2.cds_start)
    end_test   = (ivc1.cds_end is None and ivc2.cds_end is None) or\
                 (ivc1.cds_end == ivc2.cds_end)
        
    return position_test & strand_test & chrom_test & start_test & end_test


#===============================================================================
# INDEX: test suites
#===============================================================================

@attr(test="unit")
@attr(speed="slow")
class test_BPlusTree(unittest.TestCase):
 
    @classmethod
    def setUpClass(cls):
        cls.cols = [3,4,5,6,8,9,12]
        cls.bedfiles = {}
        cls.bbfiles  = {}
        for col in cls.cols:
            cls.bedfiles[col] = resource_filename("plastid","test/data/annotations/100transcripts_bed%s.bed" % col) 
            cls.bbfiles[col]  = resource_filename("plastid","test/data/annotations/100transcripts_bed%s.bb" % col)
        
        cls.chrom_sizes = { }
        for line in open(resource_filename("plastid","test/data/annotations/sacCer3.sizes")):
            chrom,size = line.strip().split("\t")
            cls.chrom_sizes[chrom] = int(size)
        
        cls.chrom_id_name = { N : K for N,K in enumerate(sorted(cls.chrom_sizes.keys(),
                                                                key = lambda x: x.lower())) }
        cls.bbs = { K : BigBedReader(cls.bbfiles[K]) for K in cls.cols }
            
    def test_magic_number(self):
        for num_cols in self.bbs:
            my_reader = self.bbs[num_cols]
            my_btree  = my_reader.bplus_tree
            my_headers = my_btree.header
            
            # test magic number parsing
            self.assertEqual(my_headers["magic"],0x78CA8C91)
            
    def test_key_size(self):
        for num_cols in self.bbs:
            my_headers = self.bbs[num_cols].bplus_tree.header
            max_len = max((len(K) for K in self.chrom_sizes))
            self.assertEqual(my_headers["key_size"],max_len)
            
    def test_num_chroms(self):
        for num_cols in self.bbs:
            my_reader = self.bbs[num_cols]
            my_btree  = my_reader.bplus_tree
            my_headers = my_btree.header
            
            # assert number of columns is correct
            self.assertEqual(my_headers["num_chroms"],len(self.chrom_sizes))
    
    def test_chrom_sizes(self):
        for num_cols in self.bbs:
            my_btree = self.bbs[num_cols].bplus_tree
            for k,v in self.chrom_sizes.items():
                self.assertEqual(v,my_btree.chrom_sizes[k])
    
    def test_chrom_id_name(self):
        for num_cols in self.bbs:
            my_btree = self.bbs[num_cols].bplus_tree
            for my_id, my_name in self.chrom_id_name.items():
                self.assertEqual(my_name,my_btree.chrom_id_name[my_id])
    
    def test_chrom_name_id(self):
        for num_cols in self.bbs:
            my_btree = self.bbs[num_cols].bplus_tree
            for my_id, my_name in self.chrom_id_name.items():
                self.assertEqual(my_id,my_btree.chrom_name_id[my_name])

@attr(test="unit")
class test_RTree(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.bbfile = resource_filename("plastid","test/data/annotations/100transcripts_bed12.bb")
        
        cls.chrom_sizes = { }
        for line in open(resource_filename("plastid","test/data/annotations/sacCer3.sizes")):
            chrom,size = line.strip().split("\t")
            cls.chrom_sizes[chrom] = int(size)
        
        cls.chrom_id_name = { N : K for N,K in enumerate(sorted(cls.chrom_sizes.keys(),
                                                                key = lambda x: x.lower())) }
        
        cls.bb = BigBedReader(cls.bbfile)
        
        cls.flybbfile = resource_filename("plastid","test/data/annotations/dmel-all-no-analysis-r5.54.bb")
        
    def test_magic_number(self):
        self.assertEqual(self.bb.r_tree.header["magic"],0x2468ACE0)
    
    def test_node_overlaps_roi(self):
        # overlap tests. Nodes represented as tuples of 
        #     (start_chrom,end_chrom,start_base,end_bse)
        #
        # Intervals represented as tuples of (start_chrom,start_base,end_base)
        
        # dict of node -> regions node should overlap
        overlap_tests = { 
                          (0,5,500,10000) : [(0,1,1000), # on first chromosome
                                             (0,0,1000),
                                             (0,100,505),
                                             (0,100,20000),
                                             (3,20000,20000),
                                             (1,1,1000), # in middle chromosomes
                                             (3,1,1000),
                                             (3,1,20000),
                                             (3,20000,20000), # in middle chromosome, but bounds outside last chrom bounds
                                             (5,1,5000), # on last chromosome, inside bounds
                                             (5,9500,10000), # on last chromosome, at edge of bounds
                                             (5,9500,10500), # on last chromosome, crossing bounds
                                          ],
                           (0,0,5000,10000) : [(0,5000,5001),
                                               (0,9999,10000),
                                               (0,6000,7000),
                                               (0,2000,6000),
                                               (0,9000,11000),
                                               (0,1000,11000),
                                               (0,5000,10000)
                                              ]
                         }
        non_overlap_tests = {
                          (0,0,5000,10000) : [(0,10000,20000),
                                              (0,10001,20000),
                                              (0,2000, 2500),
                                              (5,6000,8000),
                                              (5,6000,80000),
                                              (5,2000,6000),
                                              (5,2000,12000),
                                              ],
                         }
        
        for node_tup, node_tests in overlap_tests.items():
            node = dict(zip(["start_chrom_id","end_chrom_id","start_base","end_base"],
                            node_tup))
            for test in node_tests:
                self.assertTrue(RTree.node_overlaps_roi(node,*test),
                                "Failed positive overlap test: %s,%s" % (str(node_tup),str(test)))

        for node_tup, node_tests in non_overlap_tests.items():
            node = dict(zip(["start_chrom_id","end_chrom_id","start_base","end_base"],
                            node_tup))
            for test in node_tests:
                self.assertFalse(RTree.node_overlaps_roi(node,*test),
                                "Failed negative overlap test: %s,%s" % (str(node_tup),str(test)))
    
    def test_find_leaves_no_memorize(self):
        
        flybb = BigBedReader(self.flybbfile)
        self.assertEqual(flybb.r_tree.header["num_nodes"],73)
        self.assertEqual(len(flybb.r_tree.leaf_data),0)
        flybb.r_tree._find_leaves()
        self.assertEqual(flybb.r_tree.header["num_nodes"],len(flybb.r_tree.leaf_data))
        
        bb = BigBedReader(self.bbfile,memorize_r_tree=False)
        num_leaves = bb.r_tree.header["num_nodes"]
        self.assertEqual(num_leaves,17)
        
        # No leaves before we look
        self.assertEqual(len(bb.r_tree.leaf_data),0)
        
        bb.r_tree._find_leaves()
        
        # Make sure we find all
        self.assertEqual(len(bb.r_tree.leaf_data),num_leaves)
        
        # Make sure we didn't populate other data
        self.assertEqual(len(bb.r_tree.node_parent_offsets),0)
        self.assertEqual(len(bb.r_tree.node_child_offsets),0)
        self.assertEqual(len(bb.r_tree.node_genome_boundaries),0)
    
    def test_find_leaves_memorize(self):

        bb = BigBedReader(self.bbfile,memorize_r_tree=True)
        num_leaves = bb.r_tree.header["num_nodes"]
        
        # No leaves before we look
        self.assertEqual(len(bb.r_tree.leaf_data),0)
        bb.r_tree._find_leaves()
        
        # Make sure we find all
        self.assertEqual(len(bb.r_tree.leaf_data),num_leaves)
        
        # Make sure we populated other data
        self.assertEqual(len(bb.r_tree.node_parent_offsets),num_leaves)
        self.assertEqual(len(bb.r_tree.node_genome_boundaries),num_leaves)
        
        #self.assertEqual(len(bb.r_tree.node_child_offsets),num_nodes - num_leaves)
    
    def test_getitem(self):
        flybb = BigBedReader(self.flybbfile)
        # give an ROI and return  List of tuples of *(file_offset, bytes_to_read)*
        # hard answers obtained by manual checking in expt 258
        big_roi = GenomicSegment("2L",0,23011547,"+")
        expected = [(850, 14694),
                    (15544, 14890),
                    (30434, 15558),
                    (45992, 14665),
                    (60657, 15445),
                    (76102, 14645),
                    (90747, 15934),
                    (106681, 14673),
                    (121354, 15031),
                    (136385, 14654),
                    (151039, 14742),
                    (165781, 13265),
                    (179046, 4028)]
        self.assertListEqual(flybb.r_tree[big_roi],expected)
        
        small_roi = GenomicSegment("2L",0,10000,"+")
        expected = [(850, 14694)]
        self.assertListEqual(flybb.r_tree[small_roi],expected)
    
    def test_iter(self):
        # should fetch a list of data pointers of leaves
        # we should have 17 of them
        # and they should be sorted by order
        leaves = list(self.bb.r_tree)
        self.assertEqual(len(leaves),self.bb.r_tree.header["num_nodes"])
        for items in leaves:
            self.assertEqual(len(items),2)

        #for the fly, we should have 73
        flybb = BigBedReader(self.flybbfile)
        leaves = list(flybb.r_tree)
        self.assertEqual(len(leaves),flybb.r_tree.header["num_nodes"])
        for items in leaves:
            self.assertEqual(len(items),2)

    
@attr(test="unit")
class test_BigBedReader(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.cols = [3,4,5,6,8,9,12]
        cls.bedfiles = {}
        cls.bbfiles  = {}
        for col in cls.cols:
            cls.bedfiles[col] = resource_filename("plastid","test/data/annotations/100transcripts_bed%s.bed" % col) 
            cls.bbfiles[col]  = resource_filename("plastid","test/data/annotations/100transcripts_bed%s.bb" % col)
        
        cls.chrom_sizes = { }
        for line in open(resource_filename("plastid","test/data/annotations/sacCer3.sizes")):
            chrom,size = line.strip().split("\t")
            cls.chrom_sizes[chrom] = int(size)
        
        cls.bbs = { K : BigBedReader(cls.bbfiles[K],return_type=Transcript) for K in cls.cols }

        # comparisons against genome hash
        cls.binsize = 10000
        transcripts    = list(BED_to_Transcripts(open(cls.bedfiles[12])))

        cls.tx_dict     = {}
        cls.cds_dict    = {}
        cls.as_cds_dict = {}
        for tx in transcripts:
            txid = tx.get_name()
            cls.tx_dict[txid]     = tx
            cds_ivc = tx.get_cds()
            cds_ivc.attr["ID"] = txid
            if cds_ivc.length > 0:
                cls.cds_dict[txid]    = tx.get_cds()
                cls.as_cds_dict[txid] = tx.get_cds().get_antisense()
                cls.as_cds_dict[txid].attr["ID"] = txid

        cls.tx_hash     = GenomeHash(cls.tx_dict,do_copy=False,binsize=cls.binsize)
        cls.cds_hash    = GenomeHash(cls.cds_dict,do_copy=False,binsize=cls.binsize)
        cls.as_cds_hash = GenomeHash(cls.as_cds_dict,do_copy=False,binsize=cls.binsize)
        
        cls.shuffled_indices = list(range(len(transcripts)))
        shuffle(cls.shuffled_indices)

        cls.flybbfile = resource_filename("plastid","test/data/annotations/dmel-all-no-analysis-r5.54.bb")
        cls.flybedfile = resource_filename("plastid","test/data/annotations/dmel-all-no-analysis-r5.54.bed")
        
        # BigBed files with and without extra columns, with and without autoSql descriptions
        cls.bb_bonuscols = { "bb4as"     : resource_filename("plastid","test/data/annotations/100transcripts_bed4plus_bonus_as.bb"),
                             "bb12as"    : resource_filename("plastid","test/data/annotations/100transcripts_bed12plus_bonus_as.bb"),
                             "bb4no_as"  : resource_filename("plastid","test/data/annotations/100transcripts_bed4plus_bonus_no_as.bb"),
                             "bb12no_as" : resource_filename("plastid","test/data/annotations/100transcripts_bed12plus_bonus_no_as.bb"),
                           }
        cls.bonus_col_file = resource_filename("plastid","test/data/annotations/bonus_bed_columns.txt")
        
    def test_parse_header(self):
        for col,my_reader in self.bbs.items():
            my_headers = my_reader.header
            
            # test magic number parsing
            self.assertEqual(my_headers["magic"],0x8789F2EB)
            
            # assert number of columns is correct
            self.assertEqual(my_headers["field_count"],col)
            self.assertEqual(my_headers["bed_field_count"],col)

    def test_count_records(self):
        for _,my_reader in self.bbs.items():
            # make sure we have all records
            self.assertEqual(my_reader.num_records,100)
    
    def test_num_chroms(self):
        for _,my_reader in self.bbs.items():
            self.assertEqual(my_reader.num_chroms,17)
            
    def test_chrom_sizes(self):
        for _,my_reader in self.bbs.items():
            for k,v in self.chrom_sizes.items():
                self.assertEqual(my_reader.chrom_sizes[k],v)   
                
    def test_parse_bplus_tree_header(self):
        for _,my_reader in self.bbs.items():
            my_headers = my_reader.bplus_tree_header
            
            # test magic number parsing
            self.assertEqual(my_headers["magic"],0x78CA8C91)
            
            # make sure we have all 17 contigs
            self.assertEqual(my_headers["num_chroms"],17)
            
            # test chrom sizes
            for k,v in self.chrom_sizes.items():
                self.assertEqual(my_reader.bplus_tree.chrom_sizes[k],v)
    
    def test_iter_same_as_bed_reader(self):
        # implicitly tests iterate_over_chunk over all bed files, too
        for col in self.cols:
            bigbed = self.bbs[col]
            bed    = BED_to_Transcripts(open(self.bedfiles[col]))
               
            for n, (tx1, tx2) in enumerate(zip(bed,bigbed)):
                self.assertTrue(transcript_identical(tx1,tx2))
               
            self.assertEqual(n,100-1)
         
        flybb  = BigBedReader(self.flybbfile,return_type=Transcript)
        flybed = BED_to_Transcripts(open(self.flybedfile))
        for n, (tx1,tx2) in enumerate(zip(flybb,flybed)):
            self.assertTrue(transcript_identical(tx1,tx2))
         
        self.assertEqual(n,32682-1)
  
    def test_iterate_over_chunk(self):
        # we will have to hard-code the answers to this
        flybb = BigBedReader(self.flybbfile)
        fetched = list(flybb._iterate_over_chunk(850,14694))
        self.assertEqual(len(fetched),512)
        for ivc in fetched:
            self.assertEqual(ivc.spanning_segment.chrom,"2L")
            self.assertGreaterEqual(ivc.spanning_segment.start,7528)
            self.assertLessEqual(ivc.spanning_segment.end,1858550)
        
    def test_getitem_stranded(self):
        """Test fetching of overlapping features, minding strand
        
        1.  Make sure each feature can fetch its own subregion from its own neighborhood
        
        2.  Make sure each feature cannot fetch its own antisense subregion
        
        3.  Make sure each features fetches exactly the same features as a GenomeHash
        """             
        # make sure we can fetch each transcript's own CDS
        bb = self.bbs[12]
        u = 0
        for txid,cds in list(self.cds_dict.items())[:100]:
            gh_ol_features = self.tx_hash.get_overlapping_features(cds,stranded=True)
            bb_ol_features = bb[cds]
            self.assertIn(txid,(X.get_name() for X in gh_ol_features),
                          msg="%s failed to fetch its own CDS on correct strand" % txid)
            
            # make sure bb fetch matches GenomeHash fetch
            self.assertSetEqual(set([str(X) for X in gh_ol_features]),
                                set([str(X) for X in bb_ol_features]))
        
            u += 1
            
        self.assertGreater(u,0)

        # make sure we don't fetch each transcript's own antisense CDS
        # on opposite strand
        for txid,cds in list(self.as_cds_dict.items())[:100]:
            gh_ol_features = self.tx_hash.get_overlapping_features(cds,stranded=True)
            bb_ol_features = bb[cds]
            self.assertNotIn(txid,(X.get_name() for X in gh_ol_features),
                             msg="%s fetched its own name on wrong strand!" % txid)       
            self.assertSetEqual(set([str(X) for X in gh_ol_features]),
                                set([str(X) for X in bb_ol_features]))

    def test_getitem_unstranded(self):
        """Test fetching of overlapping features, disregarding strand
        
        1.  Make sure each feature can fetch its own subregion from its own neighborhood
        
        2.  Make sure each feature can fetch its own antisense subregion
        
        3.  Make sure each features fetches exactly the same features as a GenomeHash
        """               
        # make sure we can fetch each transcript's from its own CDS on same strand
        bb = self.bbs[12]
        u = 0
        for txid,cds in list(self.cds_dict.items())[:100]:
            gh_ol_features = self.tx_hash.get_overlapping_features(cds,stranded=False)
            bb_ol_features = bb.__getitem__(cds,stranded=False)
            self.assertIn(txid,(X.get_name() for X in gh_ol_features),
                          msg="%s failed to fetch its own CDS on same strand" % txid)
            
            # make sure bb fetch matches GenomeHash fetch
            self.assertSetEqual(set([str(X)+X.get_name() for X in gh_ol_features]),
                                set([str(X)+X.get_name() for X in bb_ol_features]))
        
            u += 1
            
        self.assertGreater(u,0)

        # make sure we can fetch each transcript's from its own antisense CDS
        # on opposite strand
        for txid,cds in list(self.as_cds_dict.items())[:100]:
            gh_ol_features = self.tx_hash.get_overlapping_features(cds,stranded=False)
            bb_ol_features = bb.__getitem__(cds,stranded=False)
            self.assertIn(txid,(X.get_name() for X in gh_ol_features),
                          msg="%s failed to fetched its own name on opposite strand!" % txid)   
            s1 = set([str(X)+X.get_name() for X in gh_ol_features])
            s2 = set([str(X)+X.get_name() for X in bb_ol_features])
            self.assertSetEqual(s1,s2,
                                msg="%s failure:\n    Only in first set: %s\n    Only in second set: %s" % (txid,
                                                                                                   s1-s2,
                                                                                                   s2-s1))
    def test_return_type(self):
        bb = self.bbs[12]
        i = iter(bb)
        for _ in range(5):
            self.assertTrue(isinstance(next(i),Transcript))
        ivcbb = BigBedReader(self.bbfiles[12],return_type=SegmentChain)
        i = iter(ivcbb)
        for _ in range(5):
            self.assertTrue(isinstance(next(i),SegmentChain))

    def test_get_autosql_str(self):
        for k in (4,12): 
            bbplus_as = BigBedReader(self.bb_bonuscols["bb%sas" % k])
            expected_as = open(resource_filename("plastid","test/data/annotations/bed%s_bonus_bed_columns.as" % k)).read()
            self.assertEqual(bbplus_as._get_autosql_str(),expected_as)

    def test_get_no_autosql_str(self):
        for k in (4,12):
            bbplus_noas = BigBedReader(self.bb_bonuscols["bb%sno_as" % k])
            self.assertEqual(bbplus_noas._get_autosql_str(),"")
    
    def test_custom_columns_names_with_autosql(self):
        expected = OrderedDict([("my_floats","some float values"),
                                ("my_sets","some set options"),
                                ("my_ints","signed integer values"),
                                ("my_strs","str representation of transcripts"),
                                ("my_colors","r,g,b colors"),
                                ])
        for k in (4,12): 
            fn = "bb%sas" % k
            bb = BigBedReader(self.bb_bonuscols[fn])
            self.assertEqual(bb.custom_fields,expected)
     
    def test_custom_columns_names_without_autosql(self):
        expected = OrderedDict([("custom_0","no description"),
                                ("custom_1","no description"),
                                ("custom_2","no description"),
                                ("custom_3","no description"),
                                ("custom_4","no description"),
                                ])
        for k in (4,12): 
            fn = "bb%sno_as" % k
            bb = BigBedReader(self.bb_bonuscols[fn])
            self.assertEqual(bb.custom_fields,expected)
    
    def test_custom_columns_retval_type_with_autosql(self):
        values = { "my_floats" : [],
                   "my_sets"   : [],
                   "my_ints"   : [],
                   "my_strs"   : [],
                   "my_colors" : [],
                  }
        bfile = open(self.bonus_col_file)
        for line in bfile:
            items = line.strip("\n").split("\t")
            values["my_floats"].append(float(items[0]))
            if items[1] == "":
                values["my_sets"].append(set())
            else:
                values["my_sets"].append(set([X.strip() for X in items[1].split(",")]))
            values["my_ints"].append(int(items[2]))
            values["my_strs"].append(items[3])
            values["my_colors"].append(tuple([int(X) for X in items[4].split(",")]))

        bfile.close()
        for k in (4,12):
            fn = "bb%sas" % k
            # ignore a Warning caused by trying to turn the BED color field
            # to an int- this has to deal with the fact that BedToBigBed wants
            # field 9 (itemRgb, typically uint[3]) to be `reserved uint;` in 
            # autoSql declarations
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                bb = BigBedReader(self.bb_bonuscols[fn])
                for n, item in enumerate(bb):
                    for key in values:
                        self.assertEqual(values[key][n],item.attr[key])
    
    def test_custom_columns_retval_type_without_autosql(self):
        values = { "custom_%s" % X : copy.deepcopy([])for X in range(5) }
        bfile = open(self.bonus_col_file)
        for line in bfile:
            items = line.strip("\n").split("\t")
            values["custom_0"].append(items[0])
            values["custom_1"].append(items[1])
            values["custom_2"].append(items[2])
            values["custom_3"].append(items[3])
            values["custom_4"].append(items[4])

        bfile.close()
        for k in (4,12):
            fn = "bb%sno_as" % k
            bb = BigBedReader(self.bb_bonuscols[fn])
            for n, item in enumerate(bb):
                for key in values:
                    self.assertEqual(values[key][n],item.attr[key])
