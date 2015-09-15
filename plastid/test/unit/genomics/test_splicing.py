#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.genomics.splicing"""

from plastid.genomics.roitools import SegmentChain
from plastid.genomics.splicing import get_junction_tuple#, get_junction_name, get_junction_seq
from nose.tools import assert_equal
from nose.plugins.attrib import attr

#===============================================================================
# INDEX: ivc -> tuple, string, et c conversion
#===============================================================================

@attr(test="unit")
def test_get_junction_tuple():
    for txid in known_juncs:
        ivcs   = known_juncs.get(txid,[])
        tuples = known_juncs_tuples.get(txid,[])
        assert_equal(len(ivcs),len(tuples))
        for ivc, tup in  zip(ivcs,tuples):
            yield check_get_junction_tuple, ivc, tup

def check_get_junction_tuple(ivc,expected_tup):
    """Check output of :py:func:`splice_to_tuple`
    
    Parameters
    ----------
    ivc : |SegmentChain|
         A two-exon fragment representing a splice junction

    expected_tup : tuple
        (chromosome name, half-open end of fiveprime exon, first position of threeprime exon, strand)
    """
    found_tup = get_junction_tuple(ivc)
    assert_equal(found_tup,expected_tup,"Expected %s, got %s." % (expected_tup,found_tup))


# def test_get_junction_name():
#     for txid in known_juncs:
#         ivcs   = known_juncs.get(txid,[])
#         names  = known_junc_names.get(txid,[])
#         assert_equal(len(ivcs),len(names))
#         for ivc, tup in  zip(ivcs,names):
#             yield check_get_junction_name, ivc, tup
#             
# def check_get_junction_name(ivc,expected_name,flank):
#     found_name = get_junction_name(ivc)
#     assert_equal(found_name,expected_name,"Expected %s, got %s." % (expected_name,found_name))
# 
# def test_get_junction_seq():
#     assert False
# 
# def check_get_junction_seq():
#     assert False

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

known_juncs_tuples = {
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