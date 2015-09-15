#!/usr/bin/env python
"""Unit tests for functions in :py:mod:`plastid.bin.slidejuncs`
"""

from nose.plugins.attrib import attr
from nose.tools import assert_equal, assert_set_equal, assert_greater# , assert_true, assert_greater_equal, assert_raises,
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from plastid.bin.slidejuncs import find_match_range, \
                                       find_known_in_range, \
                                       find_canonicals_in_range, \
                                       covered_by_repetitive
from plastid.genomics.roitools import SegmentChain, GenomicSegment
from plastid.genomics.genome_hash import GenomeHash                                       


#===============================================================================
# INDEX: verify test data above in case it got bungled
#===============================================================================

def check_data_by_sequence(expected_ivc,test_ivc):
    """Verify that two junctions produce identical sequence
    
    Parameters
    ----------
    expected_ivc : SegmentChain
        SegmentChain corresponding to known splice junction
    
    test_ivc : SegmentChain
        SegmentChain corresponding to test splice junction that is sequence-consistant
        with ``expected_ivc``
    """
    expected_seq = expected_ivc.get_sequence(seqs)
    test_seq     = test_ivc.get_sequence(seqs)
    assert_equal(expected_seq,test_seq,"Sequences don't match in test dataset for %s and %s." % (expected_seq.get_name(),test_seq.get_name()))
    
def test_data_by_sequence():
    for txid in seqs:
        expected_ivc = known_juncs.get(txid[0],None)
        if expected_ivc is not None:
            for query_ivc in query_juncs[txid]:
                yield check_data_by_sequence, expected_ivc, query_ivc


#===============================================================================
# INDEX: testing match range
#===============================================================================

def check_find_match_range(txid,max_slide):
    """Check output of :py:func:`find_match_range` against known data
    
    Parameters
    ----------
    txid : str
        ID of transcript to check
    """
    my_juncs = known_juncs.get(txid,[])
    for n,my_junc in enumerate(my_juncs):
        ref_minus, ref_plus = match_ranges[txid][n]
        ref_minus = max(-max_slide,  ref_minus)
        ref_plus  = min(max_slide+1, ref_plus)
        found_minus, found_plus = find_match_range(my_junc,seqs,max_slide)
        assert_equal(ref_minus,found_minus,"Minus match ranges don't match. Expected %s, got %s" % (ref_minus,found_minus))
        assert_equal(ref_plus,found_plus,  "Plus match ranges don't match. Expected %s, got %s" % (ref_plus,found_plus))

@attr(test="unit")
def test_find_match_range():
    # test case, plus strand has match range plus
    # test case, plus strand has no match range plus
    # test case, minus strand has match range plus
    # test case, minus strand has no match range plus
    for txid in seqs:
        for max_slide in (2,4,20):
            yield check_find_match_range, txid, max_slide


#===============================================================================
# INDEX: find known junctions within match range only when present
#===============================================================================

def check_find_known_in_range(query_junc):
    txid = query_junc.chrom
    minus_range,plus_range = find_match_range(query_junc,seqs,20)
    found_known = find_known_in_range(query_junc,minus_range,plus_range,all_known_juncs)
    
    expected    = set([str(X) for X in known_juncs.get(txid,[])])
    found_known = set([str(X) for X in found_known])
    assert_set_equal(expected,found_known,("Known junction not found for %s. "+\
                                           "Query: %s\n    Expected %s\n    " +\
                                           "Got %s\n    match range: (%s,%s)") % (txid,
                                                                                str(query_junc),
                                                                                expected,
                                                                                found_known,
                                                                                minus_range,
                                                                                plus_range))

@attr(test="unit")
def test_find_known_in_range():
    # test case, has known in plus range plus strand
    for txid in seqs:
        juncs = query_juncs.get(txid,[])
        for my_junc in juncs:
            yield check_find_known_in_range, my_junc

def check_find_none_when_none_known_in_range(query_junc):
    minus_range, plus_range = find_match_range(query_junc,seqs,20)
    found_known = find_known_in_range(query_junc,minus_range,plus_range,all_known_juncs)
    found_known = set([str(X) for X in found_known])
    assert_set_equal(found_known,set([]))

@attr(test="unit")
def test_find_none_when_none_known_in_range():
    for _,v in unmatched_query_juncs.items():
        for junc in v:
            yield check_find_none_when_none_known_in_range, junc


#===============================================================================
# INDEX: find canonical junctions within match range only when present
#===============================================================================

def check_find_canonicals_in_range(query_junc):
    txid = query_junc.chrom
    minus_range, plus_range = find_match_range(query_junc,seqs,20)
    found_can = find_canonicals_in_range(query_junc,minus_range,plus_range,seqs,canonicals[query_junc.strand])
    found_can = set([str(X) for X in found_can])
    expected  = set([str(X) for X in known_juncs.get(txid,[])])
    assert_set_equal(expected,found_can,("Known junction not found for %s. "+\
                                           "Query: %s\n    Expected %s\n    " +\
                                           "Got %s\n    match range: (%s,%s)") % (txid,
                                                                                str(query_junc),
                                                                                expected,
                                                                                found_can,
                                                                                minus_range,
                                                                                plus_range))

@attr(test="unit")
def test_find_canonicals_in_range():
    for txid in seqs:
        juncs = noncan_juncs.get(txid,[])
        for my_junc in juncs:
            yield check_find_canonicals_in_range, my_junc

def check_find_none_when_none_canonical_in_range(query_junc):
    minus_range, plus_range = find_match_range(query_junc,seqs,20)
    found_known = find_canonicals_in_range(query_junc,minus_range,plus_range,seqs,canonicals[query_junc.strand])
    found_known = set([str(X) for X in found_known])
    assert_set_equal(found_known,set([]))

@attr(test="unit")
def test_find_none_when_none_canonical_in_range():
    for junc in unmatched_noncan_query_juncs:
        yield check_find_none_when_none_canonical_in_range, junc
        

#===============================================================================
# INDEX: find repetitive junctions within match range only when present
#===============================================================================

@attr(test="unit")
def test_find_repetitive_in_range():
    pos = 0
    neg = 0
    # some should be positive, others negative
    for txid in known_juncs:
        expected   = txid in cross_hash_seqs
        if expected == True:
            pos += 1
        else:
            neg += 1
        my_query_juncs = query_juncs.get(txid,[])
        for query_junc in my_query_juncs:
            minus_range, plus_range = find_match_range(query_junc,seqs,20)
            yield check_find_repetitive_in_range, query_junc, minus_range, plus_range, expected
    
    # all negative
    for txid in unmatched_query_juncs:
        my_query_juncs = unmatched_query_juncs.get(txid,[])
        for query_junc in my_query_juncs:
            minus_range, plus_range = find_match_range(query_junc,seqs,20)
            yield check_find_repetitive_in_range, query_junc, minus_range, plus_range, False
    
    # make sure we found a bunch of each type
    assert_greater(pos,0)
    assert_greater(neg,0)
             
def check_find_repetitive_in_range(query_junc,minus_range,plus_range,expected):
    """Find repetitive blocks that overlap splice junctions
     
    Parameters
    ----------
    query_junc : SegmentChain
        splice junction
     
    expected : bool
        Whether or not ``query_junc`` overlaps a repetitive feature
    """
    assert_equal(covered_by_repetitive(query_junc,minus_range,plus_range,cross_hash),expected)

#===============================================================================
# INDEX: test data for tests above
#===============================================================================

canonicals = { "+" :  [("GT","AG"),
                       ("GC","AG")
                      ],
               "-" :  [("CT","AC"),
                        ("CT","GC")
                      ]
              }
"""Canonical splice junction sequences for each strand"""

seqs = {
        "YNL130C_fake" : "AAAAGGGCCTCAATACATGGTTTGAAAGGAAAGAGCGGTCATCACTTTGATACCTAAAGTCAGATTAGTAATTTAATGTTACAATAGTTAGTATTTATTTTTCTTGCCACAGAAGATTTAGGGTGCAATAAGATAAGCAACATACTTGTAAAGCTTCAAGTTTCCCAAACTACTCTGAGGAATAAAGAAT",
        "YNL130C"      : "AAAAGGGCCTCAATACATGGTTTGAAAGGAAAGAGCGGTCATCACTTTGATACCTAAAGTCAGATTAGTAATTTAATGTTACAATAGTTAGTATTTATTTTTCTTGCCACAGAAGATTTAGGGTGCAATAAGATAAGCAACATACTTGTAAAGCTTCAAGTTTCCCAAACTACTCTGAGGAATAAAGAAT",
        "YPL249C-A"    : "GGCTGGAGTCATTTGGGTGACTTTCTTACCCTTGTTCAAACCAATAGCGATACCTAGAGGAAAAAAATAGAATATCTTATGTTAGTAACAAGCACACTAGTACTCACGTAGTATACGAGCAAAATAAATAAACAAAGACTCCACGGAGTACTGCAACGCTGGGTGCAGGATATAAACGATGTTTGCATATCTCTGGAGCTTAATCTCATCCTGTTTGTATTTCTGTTGTTGGCTTAGTGGTTCAGTTCATTTTGCTTTTGCTCTACAGTTCATTAGTTCAGGTAAACATACCAGTCTTGACAGCCATTTTGTATTATCCTTGCTTCTCTATTCTG",

        'YBR215W_mRNA_0'  : 'ACACACAATGCAGAAACTGATGCGGTTAATTTTGCATTCTGAAACGATTTAAACGAACAGCATTAACCTCCACGACCATATTCAAACGATTGGAAATGGACCAAAAAGGTATGTTCTGGCTTTATTTTCAATTATCCGCCCTTCGTTAGGGGTGTTTTAGATCATTTTGATTAACAATTCTATTTATGATAGCAATTGTCCTTGATAACTCTAAAAGTGGTAGTAAACAGACAAAATCGAGTGGCAAAATGCAAACGCAAACAGATACTAATGCTGAAGTTCTAAACACCGATAACAGTATCAAAAAAGAAACAGGAAGTGATTCTGAGGATCTTTTCAATAAATTTTCAAACAAAAAAACGAACAGGAAAATACCCAACATCGCTGAAGAATTAGCGAAAAACAGAAACTACGTAAAAGGTGCGTCTCCGTCTCCCATTATAATTTCTGGTTCCTCTTCAACTTCTCCATCAGGACCCTCTTCCTCCAGTACAAACCCGATGGGGATACCAACGAATCGATTTAATAAAAATACTGTAGAACTATATCAGCATTCACCATCCCCCGTGATGACTACTAACAAGACTGATACAGAGGAAAAGAGACAAAACAATAGAAATATGGATAATAAAAACACCCCAGAAAGAGGATCTTCTTCCTTTGCTGCTAAGCAACTAAAAATATCATCCTTGTTAACTATATCCTCCAATGAAGATAGTAAAACTTTGCATATAAACGATACTAATGGCAACAAAAATAGTAATGCCGCCAGTAATAATATTCCATCGGCTTATGCAGAACTCCATACAGAGGGAAATAGCATTGAATCACTAATAAAACCCCCAAGCTCTCCCAGAAATAAATCTCTAACGCCCAAAGTTATCTTGCCAACACAGAATATGGATGGTACAATTGCAAAGGATCCTCATTTAGGTGATAATACACCAGGAATACTCATAGCAAAAACTAGTTCTCCCGTAAATTTGGATGTTGAAAGCACTGCACAATCTCTGGGAAAATTTAACAAGTCCACTAATTCTTTGAAAGCTGCCCTAACAAAAGCTCCTGCAGAAAAGGTCTCTTTAAAACGTAGCATTAGTTCAGTTACAAACAGTGATTCTAATATTAGCTCCAGCAAAAAGCCTACGTCTGAAAAAGCTAAAAAATCAAGTTCAGCATCAGCAATACTACCAAAGCCAACAACGACCAAGACATCAAAAAAAGCTGCTAGCAATAGTAGCGATTCTACCAGAAAAAAAAACGCCTCGAACAAGACTACATCAGCCATAAAGAAAGAGTCAAATGCCGGCTCCAAATTGAACACTGTTAAGAAGGAAAATAGCTCTTTATCTTCTATCAAAGCCACTGAGAAAGAAAAAGATAAAGGTGGAAATAGCACGGAAGCAAAAAACTCTACCAGCAACGTTAGAAAGGAACCAACTGCAAAATCCCCAAAAAGGTTAGTAGCCGCACCGACAGTCAGTCCACCAAAAATTCTACAAACGGCAGAAACCAAAGCAAAGGAACCCTCGATATTGATTGACGTTCCATTGTATCAAGCTGATACAAACGACTACCTAGACGAAAACGGTCAGGTTATTTTTAATTTATCAACTTTAATAAAAGAGAAGTATCACCCGAAAAGTAAAGAACTAGCACAACTTAAAGACTCTAAAAGAAATTTATTAATGCAATTATCCGATCATTCTAACGGCTCACTAGAAAAGGAAAAAGATGAAGAAGGAGATGTAATAGAACTTGATGATGATGAGGATATGGAAGAAGATGAAGGGGAAATAGATACAGAAACTAATACTGTAACTACCACCATCTCTCCGAAGAAGAAGTCGCATCCTATGAAAGGTAAAAATTTGATTGGTAAATATGACGTTGAAGATCCGTTCATTGATGATTCGGAGTTATTGTGGGAGGAACAGCGTGCTGCCACCAAGGACGGTTTTTTTGTTTACTTTGGTCCATTAATCGAAAAAGGTCACTATGCAAGTTTGGAGCGTGCAAATGGTACCATGAAAAGAGGTGGTGTTAAAAATAAATAAGCATTGTGTCAGAGATTTTCCAAAATCCCTGAACGACATAAAAAATGAAATAAAGTGTGCCTATATCTGAAAGTATATAGTCAGTTTATATGAAGATTTTTCCGAGCGAATGCTTTTC',
        'YHL001W_mRNA_0'  : 'GCGCAAATAAACCAAAAATGTCCACTGATTCTATTGTTAAAGCCTCCAACTGGAGATTAGTCGAAGTCGGCCGTGTCGTTTTGATCAAAAAGGGTCAATCCGCAGGTAAATTGGCTGCTATCGTTGAAATTATCGACCAAAAGAAGGTATGTTGAACCTAAAACCCACCGTGGACAAACTGAGGAGGAAATTGTAAGGAAGAGAAAGTCCCCGTATGTTCAGGGACCGCCAACGACTTCCTGATTCATGCCTGAATAGTGAAACATGTTAAACAATGAGTAGATGCTGAAATATGCAATGGGACAGAGGGCCGGAGATGTTTCTATTCTTTTTTTGCATAAACAACAGTGAAATTTCTACATATTTTGTTACTGGCAAACTCCTTCCCTGGCAGTTAATAACTCCCATACTTAAGAAGATAACTTCAGCACATGTAGCCATACCTCCTCAATTGCCAAATAATCTTTTCCTGTTTTCAAAGTATTTCAATTTACTAACAACTTTTATGAATAATTTTTTAAATTATTCGTTATAATTACTATAGGTTTTGATTGATGGTCCAAAAGCTGGTGTCCCACGTCAAGCCATCAACTTGGGTCAAGTTGTCTTAACTCCATTGACCTTTGCTCTACCAAGAGGTGCTAGAACTGCTACCGTTTCTAAGAAGTGGGCTGCTGCTGGTGTTTGCGAAAAGTGGGCTGCTTCATCTTGGGCTAAGAAGATTGCTCAACGTGAAAGACGTGCTGCTTTGACTGACTTTGAAAGATTCCAAGTTATGGTTTTAAGAAAGCAAAAGAGATACACTGTCAAGAAGGCTTTGGCTAAGGCTTAAAAAAAGAAGAATAATTCTAAAATCCATAGGTAAGTACTGAAAGCAATTTTGCGTTCCGTCAATGCATATTATATATATTAATCTTAACCATTTATGTAAACAACATATCATTTCATTTTGTTCTGGCCA',
        'YIL018W_mRNA_0'  : 'CCAAGAAACCACAAAGTTATTGAACAATGGGTATGTTTTTATTATATCGCATAATTATGGCAAATGTTATGAAGGATTCCTCTATGACTTAGATGTTTTGAATCGGTACATTTTATTTTCAGTATCCTTCTGCATATTACAAGGCAACATAGCAGCGGACAAGAAAGCTTTTTTGATGTTCGTCTTCGAAACGATTACTATCAGGGTCTTTTAAATGCTATATCGGGACTTATTGAAATTGACATATTACACTTATGAATGACGTTGGCTATCAAAGTATGAGAATAAGCCTTTCCTACTATTTTGTATCATGACAACAGGGTTCTGTCGCTTTGAATGCGGCCTTTTACTTTGCCATATTTTGATAATGAAAAAACTTTCGAGAAATATTTACTAACAGGATCGAGAATTTTTTTCCTCATTTAAACAGGTAGAGTTATTCGTAACCAAAGAAAGGGTGCTGGTTCTATCTTCACCTCTCACACCAGATTAAGACAAGGTGCTGCCAAGTTGAGAACTTTGGACTATGCTGAACGTCATGGTTACATCCGTGGTATCGTTAAGCAAATTGTCCACGACTCCGGTAGAGGTGCTCCATTGGCCAAGGTTGTTTTCCGTGACCCATACAAGTACAGATTACGTGAAGAAATCTTCATTGCTAACGAAGGTGTTCACACTGGTCAATTCATTTACGCTGGTAAAAAGGCTTCTTTGAACGTCGGTAACGTCTTGCCATTGGGTTCTGTCCCAGAAGGTACCATTGTCTCCAACGTTGAAGAAAAGCCAGGTGACAGAGGTGCCCTAGCCAGAGCTTCTGGTAACTACGTTATCATCATTGGTCACAACCCAGATGAAAACAAAACCAGAGTCAGATTACCATCCGGTGCCAAGAAGGTTATCTCTTCTGACGCCAGAGGTGTCATCGGTGTCATTGCCGGTGGTGGTAGAGTTGACAAACCATTGTTGAAGGCTGGTCGTGCTTTCCACAAGTACAGATTGAAGAGAAACTCTTGGCCAAAGACCCGTGGTGTTGCCATGAATCCAGTTGATCACCCTCACGGTGGTGGTAACCATCAACATATTGGTAAGGCTTCTACTATCTCTAGAGGTGCCGTTTCTGGTCAAAAGGCTGGTTTGATTGCTGCCAGAAGAACCGGTTTACTACGTGGTTCTCAAAAGACCCAAGATTAATCTTTTTAATTTTGGTTTCTTCCTTCTGTCATATTATTTTATCAATTTTCTTAAATATTATATAATTTAATCCGAAACGTTCCTTATAA',
        'YIL133C_mRNA_0'  : 'TGGTTTTTATTATTTCTAAAAAATATTTGTGATTATATTACAAAATGTAAAGTATTAAAGAAAAACTATTTTTAATTTCTTAGTAACCCAAAGCGGCTAATTGTTTAGCAACATCAGATTCAGCAGCAGTAGCATTAGCAGAGGCAACTTTCTTGGTGAAAGCTCTCTTCTTGGCATAGTATTCAGCGGATGAAACCTTTCTCTTTGCTTCCAATTTGGCAACAACATCTTCGTATTTCCAACCAACGGAAGTAGACAACTTACCCAAGGTGGTGTACTTTCTACCTGGCTTCAATCTCAAAACTCTCAATGCTTGTGGGACAACAACTCTCTTCTTCTTATCGTATGGAGGTGGGATACCTTCAAAGACCTTTAAACGTTCCAAAGCGGCCTTACCACGAGCAGTCTTGTGAGAAACCATACCACGAAGAGCTTTGTAGAAGATTCTAGATGGGGCTCTGAAATGGAATGGACCACGGGTCTTATTGAAAGCGGTAGCCTTTCTCAAGAAATCATGGTACTTTAATTTGTTTCTGAAAAATTCACCAGAAATGTTTAATTCTTCAGCTCTAACAACAACGATTTTTTGACCGTTTAATAATTGCTTAGCAACAACGGAAGCTAAACGACCTACTAAATGACCCTTACCTGGTTGAAGATTATTATTGATGAATGAATGGATGCTCATTCAAAAAAATGGATTTGTTAGTATTATCCTTTAAATAAAACCTGTTTTTCCCGTTCATAATTCCCATCTCTCTCTATCTCCATGTTTTAGCGTATACACAATCTACGAAGACCCTTCTGTATTCCACCCTTCTGCGGCTGTAATATTTTTATCATATTCCTTGGGGACTCGTAAACGCTGTGCACTATGAGCAAATCCTGCGGCTGGGAAAGGAAAACACATTCAATTATTAGTGCAAAATCTTACGTACCATCAATGACAACAACTGGTTCAACAGACATTTTCTCGATTTGTTCTTCACCTTCCGTTTTCAAAGAGA',
        'YIL156W_B_mRNA_0': 'GATCAAAACAAGAACAAGCAAAACGCAAAGAGATGACATTAGTAAGTACCCAAATGAGCTACTAACAACGCATCCGGTAATACTAACAAGAGAAATTGGTTAGGTAGGGAAACTAGTACACATATCCATAGATTTGGTCTTGGTATCAACGTGCCTCGCGGGGATTAAACGCAATACAGGTTTGACGCCTAAACTGGAGACTCTGGACAATCAGACTATGAGAAACTACATGAAGCGGTACTTGAACTTGGGAGAGTCTGTATACGACTATTCTGTTGCCACTTGTGGCTCCTCGACGTATTTCGCTAGGAAATAGGTACGTTTATATTCTATATATGTTGATAATATTTTGTATAACGGCACCCAGCCGGGGGCCAAAATTAATAACCGCAACTTAAAGAGCTATCG',
        'YKL006W_mRNA_0'  : 'AAGAGAGTCGTGAAAAATAAAATAAACAATGTCCACCGATTCTATTGTCAAGGCTTCTAACTGGAGATTAGTCGAAGTTGGCCGTGTGGTTTTGATCAAGAAGGGTCAATCCGCAGGTAAATTGGCCGCTATCGTCGAAATTATTGACCAAAAGAAGGTATGTTAAATCCGGAAAACCTATCATCGATTTGAGGAGGGAAATGAACGAGGAATGAGATTTAGGCTACGAAAAAGTGCTGGCAAATATCCAATCCGTTACCATGGATGAAGCTACTTGTAAACGATAAGTAAATGTGAAAATACACATTGCGATAGAGATACAAAAAAGTCTTTATTTCCCTGTTGGAGTGGCCTATAGTGGATTATTTATTTATTTTGCTAGCAGCTTACCTAGTTTCAGAAAGGTGATGACCCATATTTGCGCAGTTAATTCTTAAAGAGCCTGACCATGTTTCCTCAAATATTAAACCATGCCTTCCATCGTTTTGAAAAAATTTCGTTTACTAACAACTTTTATGAATAATTTTTTTAAATTATTTGTTATAATTATTATAGGTTTTGATTGACGGTCCAAAAGCTGGTGTTCCACGTCAAGCCATCAACTTGGGTCAAGTTGTCTTAACTCCATTGACCTTTGCTCTACCAAGAGGTGCTAGAACTGCTACCGTTTCTAAGAAGTGGGCTGCTGCTGCTGTCTGCGAAAAGTGGGCTGCTTCATCTTGGGCTAAGAAGATTGCTCAACGTGAAAGACGTGCTGCTTTGACTGACTTTGAAAGATTCCAAGTTATGGTTTTAAGAAAGCAAAAGAGATACACTGTCAAGAAGGCTTTGGCTAAGGCTTAAAATATGTACATGATCTTTAATTCTGATATATTTCGTATGTAATTTTATCTTTAACTGGTGATTCTTTTAATAAATAAACTACAATATTATAACTATTGAAAGCCTCTGTTA',
        'YMR194C_B_mRNA_0': 'TCTTTTTCATAAAGAAAGTGAGTGGATCTTCCCTGGAAGTACGAATAATGTACTTATGATTAGATTCTCTTCACGCACATACCTTTTTTTTTTTTTCAATTGCTAGTTCATACATACGTAGAAAGAAGGGCCAATTGGAAGCCTCAACTGTCCCCTGGTGTAAGCTTTCGCTGTTTCATCTTCAACTCAAGAAGCGAAGGTAGCGGACAGCAGGGTGACCGTGAATCTTTTCCATTGTCATTGTAAAATTTTGAGCAGCAGATGTACAGCTGATCTATAACTTTAGCGCATTTAGCATCATCGTACTGATGTGATAGAAGGCAATCTAAAATTAATGTGAAAAAACGTTAGTGAAAAGAAAAAAACGATCGAATCTTTTGAAGTTAATATTACAAACCTTGGATGGCGCAGGCCTCCTTCTGGCATGGATTAGACATATTTACCTCTTATTTGAACAAATCTCAACAACTAAAAAAACGCTCAATATCTTGGAAAAAAAGAACAGTAATATTCTGTCTCTCCACACCCTCTTCTGAACTCTCATATACACCCCATGGGTGAAAAGTAATATCCGTAAAAGTCTCGCTAACGGATGACGTTAACAACGCACGCCCGCTCCTTGCGAAATTTATTGTGGGGACGTCGTTCACGTGCAAAAGTCGCCAAATATTAGCGACGGGCAACATCGGAAAAACTCTGCGCCCGCGATTCCGATTGAGATGTGTAATT',
        'YPL249C_A_mRNA_0': 'CGATTAAGAATTACTATAATGATTTGTTACTTATTTTTTATTAACTAGCTTTGGGGGAGAGCCATGGAAAATAGCACTCGGTCTTGTGGCGGAGTTGGATTTGCTTGACTGTGGCCCCCTCACGCTGTTTAATGACGACGAGAGGCAGCAATGATGTTGTTCATTTCTTCGACCTTAGCCTTGGCTCTGGTGAAAGAACCCAATCTCTTCTTGGCGACCTTTCTGGCTCTCTTTTCACCGGAGTTTCTGATCAAATCGATCAATCTTCTTTCATATGGGGACAAACCGGCGATTTCTCTAACCAAAGATCTGACGAACTTGGTTCTGTTGGAGGCAGCACCCTTCTTGTAGGAGATCTTTGGGGCTGGAGTCATTTGGGTGACTTTCTTACCCTTGTTCAAACCAATAGCGATACCTAGAGGAAAAAAATAGAATATCTTATGTTAGTAACAAGCACACTAGTACTCACGTAGTATACGAGCAAAATAAATAAACAAAGACTCCACGGAGTACTGCAACGCTGGGTGCAGGATATAAACGATGTTTGCATATCTCTGGAGCTTAATCTCATCCTGTTTGTATTTCTGTTGTTGGCTTAGTGGTTCAGTTCATTTTGCTTTTGCTCTACAGTTCATTAGTTCAGGTAAACATACCAGTCTTGACAGCCATTTTGTATTATCCTTGCTTCTCTATTCTG'
        }

seqs = { K : SeqRecord(Seq(V,generic_dna)) for K,V in seqs.items()}
"""Sequences of contigs used in these tests"""


known_juncs = {
                 "YNL130C"   : ["YNL130C:0-53^145-180(-)",],
                 "YPL249C-A" : ["YPL249C-A:0-53^291-334(-)",],
                 
                'YBR215W_mRNA_0'  : ['YBR215W_mRNA_0:0-108^192-2175(+)'],
                'YHL001W_mRNA_0'  : ['YHL001W_mRNA_0:0-146^544-961(+)'],
                'YIL018W_mRNA_0'  : ['YIL018W_mRNA_0:0-30^430-1280(+)'],
                'YIL133C_mRNA_0'  : ['YIL133C_mRNA_0:0-648^938-1007(-)'],
                'YIL156W_B_mRNA_0': ['YIL156W_B_mRNA_0:0-41^103-408(+)'],
                'YKL006W_mRNA_0'  : ['YKL006W_mRNA_0:0-157^555-954(+)'],
                'YMR194C_B_mRNA_0': ['YMR194C_B_mRNA_0:0-325^397-729(-)'],
                'YNL130C_mRNA_0'  : ['YNL130C_mRNA_0:0-1204^1296-1382(-)'],
                'YPL249C_A_mRNA_0': ['YPL249C_A_mRNA_0:0-415^653-697(-)'],
               }
known_juncs = { K : [SegmentChain.from_str(X) for X in V] for K,V in known_juncs.items() }
"""Annotated splice junctions"""

all_known_juncs = []
for v in known_juncs.values():
    all_known_juncs.extend(v)
    
known_juncs_as_tuples = {
                 "YNL130C"   : [("YNL130C",53,145,"-"),],
                 "YPL249C-A" : [("YPL249C-A",53,291,"-"),],
                 
                'YBR215W_mRNA_0'  : [('YBR215W_mRNA_0',108,192,'+'),],
                'YHL001W_mRNA_0'  : [('YHL001W_mRNA_0',146,544,'+'),],
                'YIL018W_mRNA_0'  : [('YIL018W_mRNA_0',30,430,'+'),],
                'YIL133C_mRNA_0'  : [('YIL133C_mRNA_0',648,938,'-'),],
                'YIL156W_B_mRNA_0': [('YIL156W_B_mRNA_0',41,103,'+'),],
                'YKL006W_mRNA_0'  : [('YKL006W_mRNA_0',157,555,'+'),],
                'YMR194C_B_mRNA_0': [('YMR194C_B_mRNA_0',325,397,'-'),],
                'YNL130C_mRNA_0'  : [('YNL130C_mRNA_0',1204,1296,'-'),],
                'YPL249C_A_mRNA_0': [('YPL249C_A_mRNA_0',415,653,'-'),],
               }

query_juncs = {
                 "YNL130C"   : ["YNL130C:0-49^141-180(-)", # slid within range
                                "YNL130C:0-53^145-180(-)", # exact match
                                ],
                 "YPL249C-A" : [
                                "YPL249C-A:0-53^291-334(-)", # exact match
                                "YPL249C-A:0-49^287-334(-)", # slid in minus range
                                "YPL249C-A:0-54^292-334(-)", # slid in plus range
                                ],

                'YBR215W_mRNA_0': ['YBR215W_mRNA_0:0-106^190-2175(+)',
                                   'YBR215W_mRNA_0:0-107^191-2175(+)',
                                   'YBR215W_mRNA_0:0-108^192-2175(+)'],
                'YHL001W_mRNA_0': ['YHL001W_mRNA_0:0-144^542-961(+)',
                                   'YHL001W_mRNA_0:0-145^543-961(+)',
                                   'YHL001W_mRNA_0:0-146^544-961(+)',
                                   'YHL001W_mRNA_0:0-147^545-961(+)',
                                   'YHL001W_mRNA_0:0-148^546-961(+)'],
                'YIL018W_mRNA_0': ['YIL018W_mRNA_0:0-29^429-1280(+)',
                                   'YIL018W_mRNA_0:0-30^430-1280(+)',
                                   'YIL018W_mRNA_0:0-31^431-1280(+)',
                                   'YIL018W_mRNA_0:0-32^432-1280(+)',
                                   'YIL018W_mRNA_0:0-33^433-1280(+)'],
                'YIL133C_mRNA_0': ['YIL133C_mRNA_0:0-645^935-1007(-)',
                                   'YIL133C_mRNA_0:0-646^936-1007(-)',
                                   'YIL133C_mRNA_0:0-647^937-1007(-)',
                                   'YIL133C_mRNA_0:0-648^938-1007(-)',
                                   'YIL133C_mRNA_0:0-649^939-1007(-)'],
                'YIL156W_B_mRNA_0': ['YIL156W_B_mRNA_0:0-41^103-408(+)',
                                     'YIL156W_B_mRNA_0:0-42^104-408(+)',
                                     'YIL156W_B_mRNA_0:0-43^105-408(+)',
                                     'YIL156W_B_mRNA_0:0-44^106-408(+)'],
                'YKL006W_mRNA_0': ['YKL006W_mRNA_0:0-155^553-954(+)',
                                   'YKL006W_mRNA_0:0-156^554-954(+)',
                                   'YKL006W_mRNA_0:0-157^555-954(+)',
                                   'YKL006W_mRNA_0:0-158^556-954(+)',
                                   'YKL006W_mRNA_0:0-159^557-954(+)'],
                'YMR194C_B_mRNA_0': ['YMR194C_B_mRNA_0:0-325^397-729(-)',
                                     'YMR194C_B_mRNA_0:0-326^398-729(-)',
                                     'YMR194C_B_mRNA_0:0-327^399-729(-)'],
                'YPL249C_A_mRNA_0': ['YPL249C_A_mRNA_0:0-411^649-697(-)',
                                     'YPL249C_A_mRNA_0:0-412^650-697(-)',
                                     'YPL249C_A_mRNA_0:0-413^651-697(-)',
                                     'YPL249C_A_mRNA_0:0-414^652-697(-)',
                                     'YPL249C_A_mRNA_0:0-415^653-697(-)',
                                     'YPL249C_A_mRNA_0:0-416^654-697(-)']
               }

query_juncs = { K : [SegmentChain.from_str(X) for X in V] for K,V in query_juncs.items() }
"""Discovered splice junctions within match ranges of known/annotated splice junctions"""


noncan_juncs = {
                 "YNL130C"   : ["YNL130C:0-49^141-180(-)", # slid within range
                                ],
                 "YPL249C-A" : [
                                "YPL249C-A:0-49^287-334(-)", # slid in minus range
                                "YPL249C-A:0-54^292-334(-)", # slid in plus range
                                ],
               }
noncan_juncs = { K : [SegmentChain.from_str(X) for X in V] for K,V in noncan_juncs.items() }
"""Non-canonical splice junctions within match ranges of known/annotated splice junctions"""

match_ranges = {
                "YNL130C"   : [(-4,0),],
                "YPL249C-A" : [(-4,1),],
                'YBR215W_mRNA_0'  : [(-2, 0)],
                'YHL001W_mRNA_0'  : [(-2, 2)],
                'YIL018W_mRNA_0'  : [(-1, 3)],
                'YIL133C_mRNA_0'  : [(-3, 1)],
                'YIL156W_B_mRNA_0': [(0, 3)],
                'YKL006W_mRNA_0'  : [(-2, 2)],
                'YMR194C_B_mRNA_0': [(0, 2)],
                'YPL249C_A_mRNA_0': [(-4, 1)], 
                }
"""Match ranges for junctions in ``known_juncs`` above"""


unmatched_query_juncs = { "YNL130C" : ["YNL130C:0-23^145-180(-)",
                                       "YNL130C:0-53^165-180(-)",
                                       "YNL130C:0-70^141-180(-)",
                                       "YNL130C:0-49^121-180(-)", 
                                       "YNL130C:0-53^145-180(+)",      # exact match wrong strand
                                       "YNL130C_fake:0-53^145-180(-)", # exact match wrong chromosome
                                       ],
                         # these below are all 1 nucleotide outside match range
                        'YBR215W_mRNA_0'  : ['YBR215W_mRNA_0:0-105^189-2175(+)', 'YBR215W_mRNA_0:0-109^193-2175(+)'],
                        'YHL001W_mRNA_0'  : ['YHL001W_mRNA_0:0-143^541-961(+)', 'YHL001W_mRNA_0:0-149^547-961(+)'],
                        'YIL018W_mRNA_0'  : ['YIL018W_mRNA_0:0-28^428-1280(+)', 'YIL018W_mRNA_0:0-34^434-1280(+)'],
                        'YIL133C_mRNA_0'  : ['YIL133C_mRNA_0:0-644^934-1007(-)', 'YIL133C_mRNA_0:0-650^940-1007(-)'],
                        'YIL156W_B_mRNA_0': ['YIL156W_B_mRNA_0:0-40^102-408(+)', 'YIL156W_B_mRNA_0:0-45^107-408(+)'],
                        'YKL006W_mRNA_0'  : ['YKL006W_mRNA_0:0-154^552-954(+)', 'YKL006W_mRNA_0:0-160^558-954(+)'],
                        'YMR194C_B_mRNA_0': ['YMR194C_B_mRNA_0:0-324^396-729(-)', 'YMR194C_B_mRNA_0:0-328^400-729(-)'],
                        'YPL249C_A_mRNA_0': ['YPL249C_A_mRNA_0:0-410^648-697(-)', 'YPL249C_A_mRNA_0:0-417^655-697(-)']                         
                         }
unmatched_query_juncs = { K : [SegmentChain.from_str(X) for X in V] for K,V in unmatched_query_juncs.items() }
"""Query junctions with no known matches"""

unmatched_noncan_query_juncs = ["YNL130C:0-23^145-180(-)",
                                "YNL130C:0-53^165-180(-)",
                                "YNL130C:0-70^141-180(-)",
                                "YNL130C:0-49^121-180(-)", 
                                ]   
unmatched_noncan_query_juncs = [SegmentChain.from_str(X) for X in unmatched_noncan_query_juncs]
"""Query junctions without canonical splice junctions in the match range"""


repetitive_regions = [
    "YBR215W_mRNA_0:190-193(+)",   # threeprime splice site plus
    "YHL001W_mRNA_0:144-149(+)",   # fiveprime splice site plus
    "YIL133C_mRNA_0:935-940(-)",   # threeprime splice site minus
    "YMR194C_B_mRNA_0:325-328(-)", # fiveprime splice site minus
]
cross_hash = GenomeHash([SegmentChain(GenomicSegment.from_str(X)) for X in repetitive_regions])
cross_hash_seqs = { X.chrom for X in cross_hash.feature_dict.values() }