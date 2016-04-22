#!/usr/bin/env python
"""This script empirically determines which positions in a genome yield 
multimapping reads under a given set of alignment parameters. These positions
are saved in a `BED`_-formatted :term:`mask file`, so that they may be excluded
as from further analyses.

To identify multimapping regions, a genome sequence is diced into :term:`k-mers <k-mer>`
and the :term:`k-mers` aligned back to the genome. Positions in the genome that
give rise to :term:`k-mers <k-mer>` that align equally well to more than one
genomic location are then marked as multimapping..

`k` is specified by the user, as are the alignment parameters.


Output files
------------
The following files are made:

    OUTBASE_READLENGTH_MISMATCHES_crossmap.bed
        Final :term:`mask file` annotation, in `BED`_ format
    
    OUTBASE_READLENGTH_MISMATCHES_CHROMOSOME_kmers.fa
        :term:`K-mers <k-mer>` derived from chromosome `CHROMOSOME`. These files can
        be reused in subsequent runs allowing a different number of mismatches,
        using the ``--have_kmers`` option

where:

  - `OUTBASE` is a name meaningful to the user
  - `READLENGTH` is the :term:`k-mer` length chosen by the user
  - `MISMATCHES` is the number of mismatches permitted during alignment,
    also set by the user.

 .. note::
 
    ``crossmap`` internally uses `bowtie`_ and a bowtie index for alignments.
    Make sure you have these installed.
    
 .. note::

    For large genomes, it is highly recommended to convert the `BED`_-format output
    to a `BigBed`_, using Jim Kent's ``bedToBigBed`` utility as follows
    (from the terminal)::

        $ bowtie-inspect --summary BOWTIE_INDEX | grep Sequence |\\
                         cut -f2,3 | sed -e "s/\([^ ]\+\).*\\t/\\1\\t/"  >OUTFILE.sizes
        $ sort -k1,1 -k2,2n OUTBASE.bed > OUTBASE_sorted.bed
        $ bedToBigBed OUTBASE_sorted.bed OUTBASE.sizes OUTBASE_sorted.bb


    For small genomes (e.g. yeast, E. coli), this is unnecessary, and comes at a
    cost in speed.

    See https://github.com/ENCODE-DCC/kentUtils/tree/master/src/product/scripts
    for download & documentation of Kent utilities
"""
__author__ = "joshua"
import argparse
import sys
import os
import subprocess
import re
import inspect
import multiprocessing
import shutil
import functools
from plastid.util.io.filters import NameDateWriter, AbstractReader
from plastid.util.io.openers import get_short_name, argsopener
from plastid.genomics.roitools import SegmentChain, positionlist_to_segments, GenomicSegment
from plastid.util.scriptlib.help_formatters import format_module_docstring
from plastid.util.services.mini2to3 import xrange
from plastid.util.services.exceptions import MalformedFileError
from plastid.util.scriptlib.argparsers import SequenceParser, BaseParser

namepat = re.compile(r"(.*):([0-9]+)\(\+\)")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


BigBedMessage = """Crossmap complete and saved as 'OUTFILE.bed'.

    For large (e.g. mammalian) genomes, it is highly recommended to convert
    the BED-format output to a BigBed file, using Jim Kent's bedToBigBed
    utility as follows (from the terminal):
    
        $ bowtie-inspect --summary BOWTIE_INDEX | grep Sequence |\\
                         cut -f2,3 | sed -e "s/\\([^ ]\+\\).*\\t/\\1\\t/"  >OUTFILE.sizes
        $ sort -k1,1 -k2,2n OUTFILE.bed > OUTFILE_sorted.bed
        $ bedToBigBed OUTFILE_sorted.bed OUTFILE.sizes OUTFILE_sorted.bb
    
    
    See https://github.com/ENCODE-DCC/kentUtils/tree/master/src/product/scripts
    for download & documentation of Kent utilities
    
    
"""

def simulate_reads(seq_record,fh=sys.stdout,k=30):
    """Chops a DNA sequence into :term:`k-mers <k-mer>`, mimicking a sequencing run.
    Output is delivered in fasta format. Sequences are named for position of
    origin using 0-based indices.
    
    Parameters
    ----------
    seq_record : :py:class:`Bio.SeqRecord.SeqRecord`
        DNA sequence
    
    fh : file-like
        filehandle to write output 
        
    k : int, optional
        length of k-mers to generate (Default: `30`)
    """
    seq = str(seq_record.seq)

    for x in xrange(0,len(seq_record)-k+1):
        fh.write(">%s:%s(+)\n" % (seq_record.name,x))
        fh.write("%s\n" % seq[x:x+k])

    return None

class FastaNameReader(AbstractReader):
    """Returns names of sequences in a fasta file"""

    def filter(self,line):
        """Return next sequence name in a fasta file

        Parameters
        ----------
        line : str
            Line of text

        Returns
        -------
        str
            Name of next sequence, excluding prefix `'>'` and line terminator
        """
        if line.startswith(">"):
            return line[1:].rstrip()
        else:
            return self.__next__()

def revcomp_mask_chain(seg,k,offset=0):
    """Reverse-complement a single-interval mask, correcting for `offset`.
    
    Parameters
    ----------
    seg : |SegmentChain|
        Plus-strand mask, including `offset`

    k : int
        Length of k-mers

    offset : int, optional
        Offset from 5' end of read at which to map mask (Default: `0`)

    Returns
    -------
    |SegmentChain|
        Mask on minus strand corresponding to `seg`
    """
# Algorithm note:
#
#     Let
#         FW = plus-strand coordinate
#         RC = minus-strand coordinate
#     
#     Then
#         RC = FW + k - 1 - offset
#     
#     But we are given FW + offset, so:
#     
#         RC + offset = (FW + offset) + k - 1 - offset
#         RC = (FW + offset) + k - 1 - 2*offset   
    span = seg.spanning_segment
    new_offset = k - 1 - 2*offset
    ivminus = GenomicSegment(span.chrom,
                             span.start + new_offset,
                             span.end + new_offset,
                             "-")
    return SegmentChain(ivminus)

def fa_to_bed(toomany_fh,k,offset=0):
    """Create a `BED`_ file indicating genomic origins of reads in a `bowtie`_ ``toomany`` file
    
    Parameters
    ----------
    toomany_fh : file-like
        Open filehandle to fasta-formatted ``toomany`` file from `bowtie`_

    k : int
        Length of k-mers

    offset : int, optional
        Offset from 5' end of read at which to map read, if any (Default: `0`)

    Yields
    ------
    |SegmentChain|
        Plus-strand |SegmentChain| representing a repetitive region

    |SegmentChain|
        Minus-strand |SegmentChain| representing a repetitive region
    """
    last_chrom = None
    last_pos   = None
    start_pos  = None
    reader = FastaNameReader(toomany_fh)

    for n,read_name in enumerate(reader):
        chrom,pos = namepat.search(read_name).groups()
        pos = int(pos) + offset
        if chrom != last_chrom:
            if last_chrom is not None:
                plus_chain = SegmentChain(GenomicSegment(last_chrom,start_pos,last_pos+1,"+"))
                minus_chain = revcomp_mask_chain(plus_chain,k,offset)
                last_chrom = chrom
                start_pos  = pos
                last_pos   = pos
                yield plus_chain, minus_chain
            else:
                last_chrom = chrom
                start_pos  = pos
                last_pos   = pos
        else:
            delta = pos - last_pos
            if delta > 1:
                plus_chain = SegmentChain(GenomicSegment(chrom,start_pos,last_pos+1,"+"))
                minus_chain = revcomp_mask_chain(plus_chain,k,offset)
                last_pos   = pos
                start_pos  = pos
                yield plus_chain, minus_chain
            elif delta == 1:
                last_pos = pos
            else:
                msg = "k-mers are not sorted at read %s! Aborting." % read_name
                raise MalformedFileError(toomany_fh,msg,line_num=n)
    
    # export final feature
    plus_chain = SegmentChain(GenomicSegment(chrom,start_pos,last_pos+1,"+"))
    minus_chain = revcomp_mask_chain(plus_chain,k,offset)
    yield plus_chain, minus_chain

def chrom_worker(chrom_seq,args=None):
    name, seq_or_kmers = chrom_seq
    printer.write("Processing chromosome %s..." % name)
    base         = "%s_%s_%s_%s" % (args.outbase, args.read_length, args.mismatches, name)
    kmer_file    = "%s_kmers.fa"     % base
    toomany_file = "%s_multimap.fa"  % base
    bed_file     = "%s_crossmap.bed" % base
    
    if args.have_kmers == False:
        with open(kmer_file,"w") as kmer: # only do this step if args.have_kmers == False
            simulate_reads(seq_or_kmers,kmer,args.read_length)

        kmer.close()
    else:
        kmer_file = seq_or_kmers

    argdict = { "mismatches" : args.mismatches,
                "processors" : 1, 
                "bowtie"     : args.bowtie,
                "toomany"    : toomany_file,
                "kmers"      : kmer_file,
                "ebwt"       : args.ebwt,
                "null"       : os.devnull,
                }
    
    cmd  = "%(bowtie)s -m1 -a --best -f -v %(mismatches)s -p %(processors)s %(ebwt)s %(kmers)s --max %(toomany)s >%(null)s" % argdict

    printer.write("Aligning %s-mers for chromosome '%s' :\n\t'%s'" % (args.read_length,name,cmd))
    try:
        retcode = subprocess.call(cmd,shell=True)
        if retcode < 0 or retcode == 2:
            printer.write("Alignment for chromosome '%s' terminated with status %s" % (name,retcode))
        else:
            if os.path.exists(toomany_file):
                printer.write("Assembling multimappers from chromosome '%s' into crossmap..."% name)
                with argsopener(bed_file,args,"w") as bed_out:
                    for plus_chain, minus_chain in fa_to_bed(open(toomany_file),
                                                         args.read_length,
                                                         offset=args.offset):
                        bed_out.write(plus_chain.as_bed())
                        bed_out.write(minus_chain.as_bed())
                
                    bed_out.close()
            
            else:
                printer.write("Could not find multimapper source file '%s' ." % toomany_file)
    except OSError as e:
        printer.write("Alignment failed for chromosome '%s': %s" % (name,e))

    printer.write("Cleaning up chromosome '%s'..." % name)
    os.remove(toomany_file)
    if args.have_kmers == False and args.save_kmers == False:
        os.remove(kmer_file)

    return bed_file

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
    sp = SequenceParser()
    bp = BaseParser()
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[bp.get_parser(),sp.get_parser()])
    parser.add_argument("-k",dest="read_length",metavar="READ_LENGTH",
                        type=int,default=29,
                        help="K-mer length to generate from input file. "+
                             "(Default: 29)")
    parser.add_argument("--offset",type=int,default=14,
                        help="Offset from 5' end of plus-strand read at which to attribute score (Default: 14)")
    parser.add_argument("--mismatches",metavar="N",
                        type=int,default=0,
                        help="Number of mismatches tolerated in alignment. "+
                           "(Default: 0)")
    parser.add_argument("--bowtie",dest="bowtie",default="/usr/local/bin/bowtie",
                        type=str,
                        help="Location of bowtie binary (Default: ``/usr/local/bin/bowtie``)")
    parser.add_argument("--have_kmers",default=False,action="store_true",
                        help="If specified, use k-mer files from previous run. "+\
                             " In this case 'sequence_file' should be the value "+\
                             "'outbase' from the k-mer files you want to use.")
    parser.add_argument("--save_kmers",default=False,action="store_true",
                        help="Save k-mer files for reuse in a subsequent run.")
    parser.add_argument("-p","--processes",type=int,default=2,metavar="N",
                        help="Number of processes to use (should be <= number of chromosomes")
    parser.add_argument("ebwt",type=str,
                        help="Bowtie index of genome against which crossmap will be made. In most cases, should be generated from the same sequences that are in `sequence_file`.")
    parser.add_argument("outbase",type=str,
                        help="Basename for output files")

    args = parser.parse_args(argv)
    bp.get_base_ops_from_args(args)


    #filenames
    base         = "%s_%s_%s" % (args.outbase, args.read_length, args.mismatches)
    bed_file     = "%s_crossmap.bed" % base

    #if not os.path.exists(args.sequence_file):
    #    printer.write("Could not find source file: %s" % args.sequence_file)
    #    printer.write("Exiting.")
    #    sys.exit(1)

    if args.have_kmers == True:
        import glob
        kmer_files = glob.glob(args.sequence_file+"*kmers.fa")
        seq_pat = re.compile(r".*_([^_]*)_kmers.fa")
        seqs = { seq_pat.search(X).groups()[0] : X for X in kmer_files }
    else:
        seqs = sp.get_seqdict_from_args(args,index=True) 

    worker = functools.partial(chrom_worker,args=args)
    chroms = seqs.items()

    pool = multiprocessing.Pool(processes=args.processes)
    bed_filenames = pool.map(worker,chroms,1)
    pool.close()
    pool.join()
   
    with open(bed_file,"w") as fout:
        for f in sorted(bed_filenames):
            shutil.copyfileobj(open(f,"r"),fout)
            os.remove(f)

    fout.close()

    printer.write("Done.")
    printer.write(BigBedMessage.replace("OUTFILE",bed_file.replace(".bed","")).replace("BOWTIE_INDEX",args.ebwt))


if __name__ == "__main__":
    main()
