#!/usr/bin/env python
"""This script estimates :term:`P-site offsets <P-site offset>`, stratified by read length,
in a ribosome profiling dataset. To do so, :term:`metagene averages` are calculated for each
read length surrounding the start codon, mapping the reads to their fiveprime
ends. The start codon peak for each read length is heuristically identified
as the largest peak upstream of the start codon. The distance between that
peak and the start codon itself is taken to be the :term:`P-site offset` for that
read length.

Notes
------
Users should carefully examine output files to make sure these estimates are
reasonable, because if clear start codon peaks are not present in the data,
the algorithm described above will fail.  For this reason, in addition to the
:term:`P-site offsets <P-site offset>`, full metagene profiles are
exported as both tables and graphics.

Output files
------------
    OUTBASE_p_offsets.txt
        Tab-delimited text file with two columns. The first is read length,
        and the second the offset from the fiveprime end of that read length
        to the ribosomal P-site. This table can be supplied as the argument 
        for ``--offset`` when using ``--fiveprime_variable`` mapping in any
        of the other scripts in :obj:`plastid.bin`

    OUTBASE_p_offsets.svg
        Plot of metagene profiles for each read length, when reads are mapped
        to their 5' ends, :term:`P-site offsets <P-site offset>` are applied.

    OUTBASE_metagene_profiles.txt
        Metagene profiles, stratified by read length, before :term:`P-site offsets <P-site offset>`
        are applied.

    OUTBASE_K_rawcounts.txt
        Raw count vectors for each :term:`metagene` window specified in input ROI file,
        without P-site mapping rules applied, for reads of length `K`

    OUTBASE_K_normcounts.txt
        Normalized count vectors for each metagene window specified in input ROI file,
        without P-site mapping rules applied, for reads of length `K`

where `OUTBASE` is supplied by the user.
"""
import sys
import argparse
import inspect
import warnings

import matplotlib
matplotlib.use("Agg")
import numpy
import pandas as pd
import matplotlib.pyplot as plt

from collections import OrderedDict
from plastid.util.scriptlib.argparsers import get_genome_array_from_args,\
                                              get_alignment_file_parser,\
                                              get_plotting_parser,\
                                              get_colors_from_args,\
                                              get_figure_from_args
from plastid.genomics.roitools import SegmentChain
from plastid.util.io.openers import get_short_name, argsopener, NullWriter, opener
from plastid.util.io.filters import NameDateWriter
from plastid.util.scriptlib.help_formatters import format_module_docstring

warnings.simplefilter("once")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

disabled_args = ["normalize",
                 "big_genome",
                 "nibble",
                 "offset",
                 "fiveprime_variable",
                 "fiveprime",
                 "threeprime",
                 "center"]

def do_count(roi_table,ga,norm_start,norm_end,min_counts,min_len,max_len,printer=NullWriter()):
    """Calculate a :term:`metagene profile` for each read length in the dataset
    
    Parameters
    ----------
    roi_table : :class:`pandas.DataFrame`
        Table specifying regions of interest, generated
        by :py:func:`plastid.bin.metagene.do_generate`
    
    ga : |BAMGenomeArray|
        Count data
    
    norm_start : int
        Coordinate in window specifying normalization region start
    
    norm_end : int
        Coordinate in window specifying normalization region end
    
    min_counts : float
        Minimum number of counts in `window[norm_start:norm_end]`
        required for inclusion in metagene profile

    min_len : int
        Minimum read length to include
    
    max_len : int
        Maximum read length to include

    printer : file-like, optional
        filehandle to write logging info to (Default: :func:`~plastid.util.io.openers.NullWriter`)
               
    Returns
    -------
    dict
        Dictionary of :class:`numpy.ndarray` s of raw counts at each position (column)
        for each window (row)
    
    dict
        Dictionary of :class:`numpy.ndarray` s of normalized counts at each position (column)
        for each window (row), normalized by the total number of counts in that row
        from `norm_start` to `norm_end`
    
    :class:`pandas.DataFrame`
        Metagene profile of median normalized counts at each position across
        all windows, and the number of windows included in the calculation of each
        median, stratified by read length
    """
    window_size    = roi_table["window_size"][0]
    upstream_flank = roi_table["zero_point"][0]
    
    raw_count_dict  = OrderedDict()
    norm_count_dict = OrderedDict()
    for i in range(min_len,max_len+1):
        raw_count_dict[i] = numpy.ma.MaskedArray(numpy.tile(numpy.nan,(len(roi_table),window_size)))
        raw_count_dict[i].mask = numpy.tile(False,raw_count_dict[i].shape)
    
    for i,row in roi_table.iterrows():
        if i % 1000 == 0:
            printer.write("Counted %s ROIs..." % (i+1))
            
        roi    = SegmentChain.from_str(row["region"])
        mask   = SegmentChain.from_str(row["masked"])
        roi.add_masks(*mask)
        valid_mask = roi.get_masked_counts(ga).mask
        
        offset = int(round((row["alignment_offset"])))
        assert offset + roi.length <= window_size
        
        count_vectors = {}
        for k in raw_count_dict:
            count_vectors[k] = []

        for iv in roi:
            reads = ga.get_reads(iv)
            read_dict = {}
            for k in raw_count_dict:
                read_dict[k] = []

            for read in filter(lambda x: len(x.positions) in read_dict,reads):
                read_dict[len(read.positions)].append(read)
            
            for read_length in read_dict:
                count_vector = list(ga.map_fn(read_dict[read_length],iv)[1])
                count_vectors[read_length].extend(count_vector)
                
        for read_length in raw_count_dict:
            if roi.strand == "-":
                count_vectors[read_length] = count_vectors[read_length][::-1]

            raw_count_dict[read_length][i,offset:offset+roi.length]      = numpy.array(count_vectors[read_length])
            raw_count_dict[read_length].mask[i,offset:offset+roi.length] = valid_mask
    
    profile_table = { "x" : numpy.arange(-upstream_flank,window_size-upstream_flank) }
    
    printer.write("Counted %s ROIs total." % (i+1))
    for read_length in raw_count_dict:
        raw_count_dict[read_length] = numpy.ma.masked_invalid(raw_count_dict[read_length])
        denominator = raw_count_dict[read_length][:,norm_start:norm_end].sum(1)
        norm_count_dict[read_length] = (1.0*raw_count_dict[read_length].T / denominator).T
    
        norm_counts = numpy.ma.masked_invalid(norm_count_dict[read_length])
    
        with warnings.catch_warnings():
            # ignore numpy mean of empty slice warning, given by numpy in Python 2.7-3.4
            warnings.filterwarnings("ignore",".*mean of empty.*",RuntimeWarning)
            try:
                profile   = numpy.ma.median(norm_counts[denominator >= min_counts],axis=0)
            # in numpy under Python3.5, this is an IndexError instead of a warning
            except IndexError:
                profile = numpy.zeros_like(profile_table["x"],dtype=float)

        num_genes = ((~norm_counts.mask)[denominator >= min_counts]).sum(0) 
        
        profile_table["%s-mers" % read_length]         = profile
        profile_table["%s_regions_counted" % read_length] = num_genes
        
    profile_table = pd.DataFrame(profile_table)
    
    return raw_count_dict, norm_count_dict, profile_table


def main(argv=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
	argv : list, optional
		A list of command-line arguments, which will be processed
		as if the script were called from the command line if
		:py:func:`main` is called directrly.

        Default: `sys.argv[1:]`. The command-line arguments, if the script is
        invoked from the command line
    """
    alignment_file_parser = get_alignment_file_parser(disabled=disabled_args,
                                                      input_choices=["BAM"])
    plotting_parser = get_plotting_parser()

    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[alignment_file_parser,
                                              plotting_parser])
    
    parser.add_argument("--min_counts",type=int,default=10,metavar="N",
                         help="Minimum counts required in normalization region "+
                              "to be included in metagene average (Default: 10)")
    parser.add_argument("--norm_region",type=int,nargs=2,metavar="N",
                         default=(70,100),
                         help="Portion of ROI against which each individual profile"+
                              " will be normalized. Specify two integers, in nucleotide"+
                              " distance, from 5\' end of ROI. (Default: 70 100)")
    parser.add_argument("--require_upstream",default=False,action="store_true",
                        help="If supplied, the P-site offset is taken to be the distance "+
                             "between the largest peak upstream of the start codon and "+
                             "the start codon itself. Otherwise, the P-site offset is taken "+
                             "to be the distance between the largest peak in the entire ROI "+
                             "and the start codon."
                        )

    parser.add_argument("--default",type=int,default=13,
                        help="Default 5\' P-site offset for read lengths that are not present or evaluated in the dataset (Default: 13)")

    parser.add_argument("roi_file",type=str,
                        help="ROI file surrounding start codons, from ``metagene generate`` subprogram")
    
    parser.add_argument("outbase",type=str,help="Basename for output files")
    
    # set manual options
    args = parser.parse_args(argv)
    args.mapping = "fiveprime"
    args.offset  = 0
    args.nibble  = 0

    
    # process arguments
    min_len = args.min_length
    max_len = args.max_length
    profiles = max_len + 1 - min_len
    lengths = list(range(min_len,max_len+1))
    outbase = args.outbase
    title  = "Fiveprime read offsets by length" if args.title is None else args.title
    colors = get_colors_from_args(args,profiles)
 
    printer.write("Opening ROI file %s..." % args.roi_file)
    with opener(args.roi_file) as roi_fh:
        roi_table = pd.read_table(roi_fh,sep="\t",comment="#",index_col=None,header=0)
        roi_fh.close()
        
    printer.write("Opening count files %s..." % ",".join(args.count_files))
    ga = get_genome_array_from_args(args,printer=printer,disabled=disabled_args)

    
    # remove default size filters
    my_filters = ga._filters.keys()
    for f in my_filters:
        ga.remove_filter(f)

    # count
    count_dict, norm_count_dict, metagene_profile = do_count(roi_table,
                                                             ga,
                                                             args.norm_region[0],
                                                             args.norm_region[1],
                                                             args.min_counts,
                                                             min_len,
                                                             max_len,
                                                             printer=printer)
    
    # save counts
    profile_fn = "%s_metagene_profiles.txt" % outbase
    with argsopener(profile_fn,args,"w") as metagene_out:
        metagene_profile.to_csv(metagene_out,
                                sep="\t",
                                header=0,
                                index=False,
                                na_rep="nan",
                                columns=["x"]+["%s-mers" % X for X in lengths])
        metagene_out.close()

    printer.write("Saving raw and normalized counts...")
    for k in count_dict:
        count_fn     = "%s_%s_rawcounts.txt.gz"  % (outbase,k)
        normcount_fn = "%s_%s_normcounts.txt.gz" % (outbase,k)
        numpy.savetxt(count_fn,count_dict[k],delimiter="\t")
        numpy.savetxt(normcount_fn,norm_count_dict[k],delimiter="\t")
        
    
    # plotting & offsets
    printer.write("Plotting and determining offsets...")
    offset_dict = OrderedDict() 

    # Determine scaling factor for plotting metagene profiles
    max_y = numpy.nan 
    with warnings.catch_warnings():
        # ignore warnings for slices that contain only NaNs
        warnings.simplefilter("ignore",category=RuntimeWarning)
        for k in lengths:
            max_y = numpy.nanmax([max_y,
                                  numpy.nanmax(metagene_profile["%s-mers"% k].values)])

    if numpy.isnan(max_y):
        max_y = 1.0


    # parse arguments & set styles
    try:
        import matplotlib.style
        if getattr(args,"stylesheet",None) is not None:
            matplotlib.style.use(args.stylesheet)
    except ImportError:
        pass

    mplrc = matplotlib.rcParams
    plt_incr  = 1.2

    # use this figsize if not specified on command line
    figheight = 1.0 + 0.25*(profiles-1) + 0.75*(profiles)
    default_figsize = (7.5,figheight)

    fig = get_figure_from_args(args,figsize=default_figsize)

    ax = plt.gca()
    plt.title(title)
    plt.xlabel("Distance from CDS start, (nt; 5' end mapping)")
    plt.ylabel("Median normalized read density (au)")
    plt.axvline(0.0,color=mplrc["axes.edgecolor"],dashes=[3,2])

    x = metagene_profile["x"].values
    xmin = x.min()
    xmax = x.max()
    mask = numpy.tile(True,len(x)) if args.require_upstream == False else (x <= 0)

    for n,k in enumerate(lengths):
        color = colors[n]
        baseline = plt_incr*n
        y = metagene_profile["%s-mers" % k].values
        ymask = y[mask]

        if numpy.isnan(y).all():
            plot_y = numpy.zeros_like(x)
        else:
            plot_y = y / max_y
 
        # plot metagene profiles on common scale, offset by baseline from bottom to top
        ax.plot(x,baseline + plot_y,color=color)
        ax.text(xmin,baseline,"%s-mers" % k,
                ha="left",
                va="bottom",
                color=color,
                transform=matplotlib.transforms.offset_copy(ax.transData,fig,
                                                            x=6.0,y=3.0,units="points"))

        ymax = baseline + plot_y.max()

        if mask.sum() == numpy.isnan(ymask).sum() or ymask.sum() == 0:
            offset = args.default
            usedefault = True
        else:
            offset = -x[ymask.argmax()]
            usedefault = False

        offset_dict[k] = offset
        if usedefault == False:
            yadj = ymax - 0.2 * plt_incr

            ax.plot([-offset,0],[yadj,yadj],color=color,dashes=[3,2])
            ax.text(-offset / 2.0,
                     yadj,
                     "%s nt" % (offset),
                     color=color,
                     ha="center",
                     va="bottom",
                     transform=matplotlib.transforms.offset_copy(ax.transData,fig,
                                                                 x=0.0,y=3.0,units="points")
                    )   

    plt.xlim(xmin,xmax)
    plt.ylim(-0.1,plt_incr+baseline)
    ax.yaxis.set_ticks([])

    # save data as p-site offset table
    fn = "%s_p_offsets.txt" % outbase
    fout = argsopener(fn,args)
    printer.write("Writing offset table to %s ..." % fn)
    fout.write("length\tp_offset\n")
    for k in offset_dict:
        fout.write("%s\t%s\n" % (k,offset_dict[k]))
    
    fout.write("default\t%s" % args.default)
    
    fout.close()

    # save plot
    plot_fn ="%s_p_offsets.%s" % (outbase,args.figformat) 
    printer.write("Saving plot to %s ..." % plot_fn)
    plt.savefig(plot_fn,dpi=args.dpi,bbox_inches="tight")

    printer.write("Done.")

if __name__ == "__main__":
    main()
