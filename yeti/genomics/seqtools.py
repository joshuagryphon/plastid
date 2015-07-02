#!/usr/bin/env python
"""Utilities for mutating and searching nucleic acid sequences"""
import random, re

IUPAC_TABLE = { "A" : "A",
                "C" : "C",
                "T" : ("T","U"),
                "U" : ("T","U"),
                "G" : "G",
                "N" : ("A","C","T","G","U"),
                "R" : ("A","G"),     # puRines
                "Y" : ("C","T","U"), # pYrimidines
                "S" : ("G","C"),     # Strong binding
                "W" : ("A","T","U"), # Weak binding
                "K" : ("G","T","U"),
                "M" : ("A","C"),
                "B" : ("C","G","T","U"),
                "D" : ("A","G","T","U"),
                "H" : ("A","C","T","U"),
                "V" : ("A","C","G"),
               }
"""Dictionary mapping IUPAC nucleotide symbols to tuples of nucleotides
they represent (e.g. R -> (A, G) )
"""

def seq_to_regex(inp,flags=0):
    """Convert a nucleotide sequence of IUPAC nucleotide characters as a regular expression.
    Ambiguous IUPAC characters are converted to groups (e.g. `'Y'` to `'[CTU]'`),
    and T and U are considered equivalent.
    
    Parameters
    ----------
    inp : str
        Nucleotide sequence using IUPAC nucleotide codes
    
    flags : int, optional
        Flags to pass to :py:func:`re.compile` (Default: 0 / no flags)
        
    Examples
    --------
    Convert a sequence to a regex::
    
        >>> seq_to_regex("CARYYA").pattern
        'CA[AG][CTU][CTU]A'
    
    
    Returns
    -------
    :py:class:`re.RegexObject`
        Regular expression pattern corresponding to IUPAC sequence in `inp`
    """
    out = []
    for ch in inp:
        if len(IUPAC_TABLE.get(ch,ch)) == 1:
            out.append(ch)
        else:
            out.append("["+"".join(IUPAC_TABLE.get(ch,ch))+"]")
    
    return re.compile("".join(out),flags=flags)

def mutate_seqs(seqs,nucleotides="NACTG",mutations=1):
    """Generate all sequences within `mutations` distance from a reference sequence
    
    Parameters
    ----------
    seqs : str or list of str
        Single reference sequence (a string) or a group of strings
                        
    nucleotides : list of char, optional
        Permitted nucleotide substitutions (Default: `'NACTG'`)
    
    mutations : int, optional
        Number of substitutions to make (Default: `1`)
    
    
    Returns
    -------
    set
        all sequences within `mutations` substitutions from the sequence(s)
        specified in `seqs`
    """
    if isinstance(seqs,str):
        seqs = [seqs]
    if mutations == 0:
        return set(seqs)
    else:
        seqsout = []
        for seq in seqs:
            for nuc in nucleotides:
                for i in range(len(seq)):
                    newseq = list(seq)[:]
                    newseq[i] = nuc
                    seqsout.append("".join(newseq))
        seqsout.extend(mutate_seqs(seqsout,nucleotides=nucleotides,mutations=mutations-1))
        return set(seqsout) | set(seqs)


def random_seq(size,nucleotides="ACTG"):
    """Generate a random nucleotide sequence of length `size` and composition `nucleotides`
    
    Parameters
    ----------
    size : int
        length of desired sequence
        
    nucleotides : str, optional
        string of nucleotides to use in sequence, in desired base composition
        (i.e. need not be unique; can supply `'AATCG'` to increase `'A'` bias.
        Default: `'ACTG'`)
        
    Returns
    -------
    str : randomized DNA sequence
    """
    seq = "".join([nucleotides[random.randrange(0,len(nucleotides))] for _ in range(0,size)])
    return seq

