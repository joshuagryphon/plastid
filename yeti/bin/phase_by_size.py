#!/usr/bin/env python
"""Estimate sub-codon phasing in a ribosome profiling dataset, stratified by read length.
NOTE: to avoid double-counting of codons, users should supply an annotation file
that includes only one transcript isoform per gene.


Output files
------------
    ${OUTBASE}_phasing.txt
        Read phasing for each read length

    ${OUTBASE}_phasing.svg
        Plot of phasing by read length

where ${OUTBASE} is supplied by the user.
"""
import sys
import matplotlib
matplotlib.use("Agg")
from yeti.util.scriptlib.argparsers import get_genome_array_from_args,\
                                                      get_transcripts_from_args,\
                                                      get_alignment_file_parser,\
                                                      get_annotation_file_parser
from yeti.util.array_table import ArrayTable
from yeti.util.io.openers import get_short_name, argsopener
from yeti.util.io.filters import NameDateWriter
from yeti.genomics.genome_array import SizeFilterFactory
from yeti.util.scriptlib.help_formatters import format_module_docstring

import numpy
import argparse
import matplotlib.pyplot as plt
import inspect

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))



def main(argv=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
	argv : list, optional
		A list of command-line arguments, which will be processed
		as if the script were called from the command line if
		:py:func:`main` is called directly.

		Default: sys.argv[1:] (actually command-line arguments)
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
        cds_ivc = roi.get_cds()
        
        # for each iv, fetch reads, sort them, and create individual count vectors
        for iv in cds_ivc:
            reads = gnd.get_reads(iv)
            for read in filter(lambda x: len(x.positions) in read_dict,reads):
                read_dict[len(read.positions)].append(read)

            for read_length in read_dict:
                count_vector = list(gnd.map_fn(read_dict[read_length],iv)[1])
                count_vectors[read_length].extend(count_vector)
        
        # add each count vector to total
        for k, vec in count_vectors.items():
            counts = numpy.array(vec)
            if roi.strand == "-":
                counts = counts[::-1]
            
            try:
                counts = counts.reshape((len(vec)/3.0,3))
                phase_sums[k] += counts[args.codon_buffer:-args.codon_buffer,:].sum(0)
            except ValueError:
                printer.write("Ignoring %s because CDS is not divisible by 3: %s." % (roi.get_name(),cds_ivc.get_length()) )

    printer.write("Counted %s ROIs total." % (n+1))
    dtmp = ArrayTable(dtmp)
    
    # total reads counted for each size    
    for k in read_lengths:
        dtmp["reads_counted"][dtmp["read_length"] == k] = phase_sums[k].sum() 
    
    # read length distribution
    dtmp["fraction_reads_counted"] = dtmp["reads_counted"].astype(float) / dtmp["reads_counted"].sum() 
    
    
    # phase vectors
    phase_vectors = { K : V.astype(float)/V.astype(float).sum() for K,V in phase_sums.items() }
    for i in range(3):
        dtmp["phase%s" % i] = numpy.zeros(len(dtmp)) 

    for k, vec in phase_vectors.items():
        for i in range(3):
            dtmp["phase%s" % i][dtmp["read_length"] == k] = vec[i]
    
    # phase table
    fn = "%s_phasing.txt" % args.outbase
    printer.write("Saving phasing table to %s ..." % fn)
    with argsopener(fn,args) as fh:
        dtmp.to_file(fh,keyorder=["read_length",
                                                   "reads_counted",
                                                   "fraction_reads_counted",
                                                   "phase0",
                                                   "phase1",
                                                   "phase2",
                                                    ],
                      formatters={numpy.float64 : '{:.6f}'.format,
                                  numpy.float32 : '{:.6f}'.format,
                                  numpy.float128 : '{:.6f}'.format,
                                  }
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
