"""Test suite for :py:mod:`plastid.genomics.seqtools`"""
import re
import sys
from nose.plugins.attrib import attr
from nose.tools import assert_equal, assert_set_equal, assert_true
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from plastid.test.ref_files import REF_FILES
from plastid.genomics.seqtools import seq_to_regex, mutate_seqs, random_seq,\
                                   _TwoBitSeqProxy, TwoBitSeqRecordAdaptor

#===============================================================================
# INDEX: unit tests
#===============================================================================

@attr(test="unit")
def test_seq_to_regex():
    for flags, inp, outp in SEQ_REGEX:
        check = seq_to_regex(inp,flags=flags)
        pattern = check.pattern
        pflags  = check.flags
        yield assert_equal, pattern, outp, "seq_to_regex(%s,%s) regex do not match. Expected %s, found %s" % (inp, flags, outp, pattern)
        yield assert_equal, flags, pflags, "seq_to_regex(%s,%s) flags do not match. Expected %s, found %s" % (inp, flags, flags, pflags)

@attr(test="unit")
def test_mutate_seqs():
    for inp, nucs, mutations, outp in MUTANTS:
        check = mutate_seqs(inp,nucleotides=nucs,mutations=mutations)
        yield assert_set_equal, check, outp, "mutate_seqs(%s,%s) did not match expectation. Expected %s, found %s" % (inp, mutations, outp, check)

def check_random_seq(length,nucleotides):
    check = random_seq(length,nucleotides=nucleotides)
    
    # check length
    assert_equal(length,len(check),"random_chr(%s,%s) length is incorrect. Expected %s, found %s" % (length,nucleotides,length,len(check)))
    
    # make sure nucleotides used are a subset of those expected
    diff = set(list(check))-set(nucleotides)
    assert_equal(diff,set(),"random_chr(%s,%s) found extra nucleotides in result: %s" % (length,nucleotides,diff))
    
@attr(test="unit")
def test_random_chr():
    for length, nucs in RANDOM_CHR:
        yield check_random_seq, length, nucs

@attr(test="unit")
def testTwoBitSeqProxyFetch():
    g1 = SeqIO.to_dict(SeqIO.parse(REF_FILES["yeast_fasta"],"fasta"))
    g2 = TwoBitSeqRecordAdaptor(REF_FILES["yeast_twobit"])
    
    seq = g2["chrV"][50:5000]
    assert_true(isinstance(seq,SeqRecord))
    assert_equal(str(g1["chrV"][50:5000].seq),str(seq.seq))

@attr(test="unit")
def testTwoBitSeqProxyLen():
    g1 = SeqIO.to_dict(SeqIO.parse(REF_FILES["yeast_fasta"],"fasta"))
    g2 = TwoBitSeqRecordAdaptor(REF_FILES["yeast_twobit"])

    for k in g1:
        assert_equal(len(g1[k]),len(g2[k]))

@attr(test="unit")
def testTwoBitSeqProxyStr():
    g1 = SeqIO.to_dict(SeqIO.parse(REF_FILES["yeast_fasta"],"fasta"))
    g2 = TwoBitSeqRecordAdaptor(REF_FILES["yeast_twobit"])

    for k in g1:
        assert_equal(str(g1[k].seq),str(g2[k]))

@attr(test="unit")
def testTwoBitSeqProxySeqProp():
    g1 = SeqIO.to_dict(SeqIO.parse(REF_FILES["yeast_fasta"],"fasta"))
    g2 = TwoBitSeqRecordAdaptor(REF_FILES["yeast_twobit"])

    for k in g1:
        assert_equal(str(g1[k].seq),str(g2[k].seq))

@attr(test="unit")
def testTwoBitSeqProxyRevComp():
    g1 = SeqIO.to_dict(SeqIO.parse(REF_FILES["yeast_fasta"],"fasta"))
    g2 = TwoBitSeqRecordAdaptor(REF_FILES["yeast_twobit"])

    for k in g1:
        assert_equal(str(g1[k].reverse_complement().seq),str(g2[k].reverse_complement().seq))

@attr(test="unit")
def testTwoBitSeqRecordAdaptor():
    g1 = SeqIO.to_dict(SeqIO.parse(REF_FILES["yeast_fasta"],"fasta"))
    g2 = TwoBitSeqRecordAdaptor(REF_FILES["yeast_twobit"])

    for k,v in g1.items():
        assert_true(isinstance(g2[k],_TwoBitSeqProxy))
   



#===============================================================================
# INDEX: test data
#===============================================================================

if sys.version_info > (3,):
    default_flag = re.UNICODE
else:
    default_flag = 0

SEQ_REGEX = [ (default_flag,"NNNN","[ACTGU][ACTGU][ACTGU][ACTGU]"),  # N->ACTGU
              (default_flag,"CCCAGA","CCCAGA"),                      # No change
              (default_flag,"TCTAGA","[TU]C[TU]AGA"),                # T->TU
              (default_flag,"UCUAGA","[TU]C[TU]AGA"),                # U->TU
              (default_flag,"ARCGA","A[AG]CGA"),                     # R->AG
              (default_flag,"YAGAC","[CTU]AGAC"),                    # Y->CTU
              (default_flag,"SAGAC","[GC]AGAC"),                     # S->GC
              (default_flag,"WAGAC","[ATU]AGAC"),                    # W->AT
              (default_flag,"KAGAC","[GTU]AGAC"),
              (default_flag,"MAGAC","[AC]AGAC"),
              (default_flag,"BAGAC","[CGTU]AGAC"),
              (default_flag,"DAGAC","[AGTU]AGAC"),
              (default_flag,"HAGAC","[ACTU]AGAC"),
              (default_flag,"VAGAC","[ACG]AGAC"),
              (default_flag,"SWANKY","[GC][ATU]A[ACTGU][GTU][CTU]"),
              
              (default_flag|re.IGNORECASE,"KAGAC","[GTU]AGAC"),
              (default_flag|re.DEBUG,"MAGAC","[AC]AGAC"),
              (default_flag|re.MULTILINE,"BAGAC","[CGTU]AGAC"),
              (default_flag|re.DEBUG|re.MULTILINE,"DAGAC","[AGTU]AGAC"),
              (default_flag|re.DOTALL|re.UNICODE,"HAGAC","[ACTU]AGAC"),
             ]
"""Test cases for :py:func:`seq_to_regex`, as tuples of
    1. Regex flags
    2. Input string
    3. Expected regex pattern
"""


MUTANTS = [("A","NACTG",1,{"A","C","T","G","N"}),
           ("C","NACTG",1,{"A","C","T","G","N"}),
           ("T","NACTG",1,{"A","C","T","G","N"}),
           ("G","NACTG",1,{"A","C","T","G","N"}),
           ("N","NACTG",1,{"A","C","T","G","N"}),
           # 2 mutations
           ("A","NACTG",2,{"A","C","T","G","N"}),
           ("AT","NACTG",1,{"AT","AA","AC","AG","AN","CT","GT","TT","NT"}),
           # 2 mutations
           ("AT","NACTG",2,{"AA","AC","AT","AG","AN",
                            "CA","CC","CT","CG","CN",
                            "GA","GC","GT","GG","GN",
                            "TA","TC","TT","TG","TN",
                            "NA","NC","NT","NG","NN"}),
           ("AAA","NACTG",1,{"AAA",
                             "AAT","AAC","AAG","AAN",
                             "ATA","ACA","AGA","ANA",
                             "TAA","CAA","GAA","NAA"}
            ),
           ("AAAA","NACTG",1,{"AAAA",
                              "AAAT","AAAC","AAAG","AAAN",
                              "AATA","AACA","AAGA","AANA",
                              "ATAA","ACAA","AGAA","ANAA",
                              "TAAA","CAAA","GAAA","NAAA",
                              }
            ),
           # restricted nucleotides
           ("AAAA","CT",1,{"AAAA",
                           "AAAT","AAAC",
                           "AATA","AACA",
                           "ATAA","ACAA",
                           "TAAA","CAAA",
                           }
            ),
           # multiple input sequences
           (["A","AAA"],"NACTG",1,{"A","C","T","G","N",
                                   "AAA",
                                   "AAT","AAC","AAG","AAN",
                                   "ATA","ACA","AGA","ANA",
                                   "TAA","CAA","GAA","NAA"}
            )   
           ]
"""Test cases for :py:func:`mutate_seqs` as tuples of:
    1. Input sequence(s)
    2. Permitted nucleotide substitutions
    3. Number of mutations
    4. Expected output
"""


RANDOM_CHR = [(1,"ACTG"),
              (100,"AC"),
              (100,"ACTG"),
              (10000,"ACTG"),
              ]
"""Test cases for :py:func:`random_seq` as tuples of:
    1. Sequence length
    2. Nucleotide composition
"""
