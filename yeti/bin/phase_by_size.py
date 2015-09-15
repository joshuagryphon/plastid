#!/usr/bin/env python
"""Estimate :term:`sub-codon phasing` in a :term:`ribosome profiling` dataset,
stratified by read length.

Because ribosomes step three nucleotides in each cycle of translation elongation,
in many :term:`ribosome profiling` datasets a triplet periodicity is observable
in the distribution of :term:`ribosome-protected footprints <footprint>`,
in which 70-90% of the reads on a codon fall within the first of the three codon
positions. This allows deduction of translation reading frames, if the reading
frame is not known *a priori.* See :cite:`Ingolia2009` for more details

Output files
------------
    OUTBASE_phasing.txt
        Read phasing for each read length

    OUTBASE_phasing.svg
        Plot of phasing by read length

where ${OUTBASE} is supplied by the user.

Notes
-----
To avoid double-counting of codons, users should supply an :term:`annotation`
file that includes only one transcript isoform per gene.

"""
import sys
import warnings
import argparse
import inspect
import warnings

import pandas as pd
import numpy
import matplotlib
matplotlib.use("Agg")
from plastid.util.scriptlib.argparsers import get_genome_array_from_args,\
                                                      get_transcripts_from_args,\
                                                      get_alignment_file_parser,\
                                                      get_annotation_file_parser
from plastid.util.io.openers import get_short_name, argsopener
from plastid.util.io.filters import NameDateWriter
from plastid.util.scriptlib.help_formatters import format_module_docstring
from plastid.util.services.exceptions import DataWarning

import matplotlib.pyplot as plt

warnings.simplefilter("once")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

# TODO: support bowtie or wiggle files

def main(argv=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
	argv : list, optional
		A list of command-line arguments, which will be processed
		as if the script were called from the command line if
		:py:func:`main` is called directly.

        Default: `sys.argv[1:]`. The command-line arguments, if the script is
        invoked from the command line
    """
    alignment_file_parser = get_alignment_file_parser(disabled=["normalize",
                                                                "big_genome",
                                                                "spliced_bowtie_files"],
                                                      input_choices=["BAM"])
    annotation_file_parser = get_annotation_file_parser()
    
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[annotation_file_parser,
                                              alignment_file_parser])
    
    parser.add_argument("--codon_buffer",type=int,default=5,
                        help="Codons before and after start codon to ignore (Default: 5)")
    
    parser.add_argument("outbase",type=str,help="Basename for output files")
    args = parser.parse_args(argv)
    gnd = get_genome_array_from_args(args,printer=printer,disabled=["normalize",
                                                                    "big_genome",
                                                                    "spliced_bowtie_files"])
    transcripts = get_transcripts_from_args(args,printer=printer)
    
    read_lengths = list(range(args.min_length,args.max_length+1))
    codon_buffer = args.codon_buffer
    
    dtmp = { "read_length"   : numpy.array(read_lengths),
             "reads_counted" : numpy.zeros_like(read_lengths,dtype=int),
            }
    
    phase_sums = {}
    for k in read_lengths:
        phase_sums[k] = numpy.zeros(3)
    
    for n, roi in enumerate(transcripts):
        if n % 1000 == 1:
            printer.write("Counted %s ROIs..." % n)
            
        read_dict     = {}
        count_vectors = {}
        for k in read_lengths:
            read_dict[k]     = []
            count_vectors[k] = []
        cds_chain = roi.get_cds()
        
        # for each seg, fetch reads, sort them, and create individual count vectors
        for seg in cds_chain:
            reads = gnd.get_reads(seg)
            for read in filter(lambda x: len(x.positions) in read_dict,reads):
                read_dict[len(read.positions)].append(read)

            # map and sort by length
            for read_length in read_dict:
                count_vector = list(gnd.map_fn(read_dict[read_length],seg)[1])
                count_vectors[read_length].extend(count_vector)
        
        # add each count vector for each length to total
        for k, vec in count_vectors.items():
            counts = numpy.array(vec)
            if roi.strand == "-":
                counts = counts[::-1]
           
            if len(counts) % 3 == 0:
                counts = counts.reshape((len(vec)/3,3))
            else:
                message = "Length of '%s' coding region (%s nt) is not divisible by 3. Ignoring last partial codon." % (len(counts),roi.get_name())
                warnings.warn(message,DataWarning)
                counts = counts.reshape(len(vec)//3,3)

            phase_sums[k] += counts[codon_buffer:-codon_buffer,:].sum(0)

    printer.write("Counted %s ROIs total." % (n+1))
    for k in dtmp:
        dtmp[k] = numpy.array(dtmp[k])
    
    # total reads counted for each size    
    for k in read_lengths:
        dtmp["reads_counted"][dtmp["read_length"] == k] = phase_sums[k].sum() 
    
    # read length distribution
    dtmp["fraction_reads_counted"] = dtmp["reads_counted"].astype(float) / dtmp["reads_counted"].sum() 
    
    # phase vectors
    phase_vectors = { K : V.astype(float)/V.astype(float).sum() for K,V in phase_sums.items() }
    for i in range(3):
        dtmp["phase%s" % i] = numpy.zeros(len(dtmp["read_length"])) 

    for k, vec in phase_vectors.items():
        for i in range(3):
            dtmp["phase%s" % i][dtmp["read_length"] == k] = vec[i]
    
    # phase table
    fn = "%s_phasing.txt" % args.outbase
    printer.write("Saving phasing table to %s ..." % fn)
    dtmp = pd.DataFrame(dtmp)
    with argsopener(fn,args) as fh:
        dtmp.to_csv(fh,columns=["read_length",
                                "reads_counted",
                                "fraction_reads_counted",
                                "phase0",
                                "phase1",
                                "phase2",
                                ],
                       float_format="%.6f",
                       na_rep="nan",
                       sep="\t",
                       index=False,
                       header=True
                      )
        fh.close()
    
    # plot
    # TODO: plot color next to line
    #       or change to bar graphs
    fn = "%s_phasing.svg" % args.outbase
    printer.write("Plotting to %s ..." % fn)
    plt.figure()
    plt.xlabel("Codon position")
    plt.ylabel("Fraction of reads")
    for i in range(len(dtmp)):
        l = dtmp["read_length"][i]
        phasing = phase_vectors[l]
        plt.plot(phasing,label="%smers" % l)
        
    plt.legend()
    plt.savefig(fn)
    

if __name__ == "__main__":
    main()
