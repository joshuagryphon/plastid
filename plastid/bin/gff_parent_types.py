#!/usr/bin/env python
"""Exports a table of parent-child feature relationships for all feature types
found in a `GFF3`_ file. Features with multiple parents are dumped into a category
called 'Multiple' and the feature types of their individual parents ignored.
"""
from plastid.readers.gff import GFF3_Reader
from plastid.util.io.filters import NameDateWriter
from plastid.util.io.openers import get_short_name, opener, argsopener
from plastid.util.scriptlib.help_formatters import format_module_docstring
from collections import Counter
import argparse
import sys
import inspect
import warnings

warnings.simplefilter("once")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

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
	parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
									 formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("--exclude",nargs="+",default=[],
					    help="Feature types to exclude from consideration")
	parser.add_argument("infile",metavar="infile.gff",type=str,
	                   help="Input GFF3 file")
	parser.add_argument("outfile",metavar="outfile.txt",type=str,
					   help="Name of output file")
	
	args = parser.parse_args(argv)
	excluded = set(args.exclude)

	fin = sys.stdin if args.infile == "-" else opener(args.infile)
	
	feature_counts        = Counter()
	features_with_parents = []
	feature_types         = {}
	name_type             = {}

	
	printer.write("Opening %s..." % args.infile)
	c = 0
	for feature in GFF3_Reader(fin,return_stopfeatures=False):
		if c % 10000 == 0:
			printer.write("Processed %s features..." % c)
		c += 1
		ftype = feature.attr["type"]
		fname = feature.get_name()
		if ftype not in excluded:
			if ftype not in feature_types:
				feature_types[ftype] = Counter()
			feature_counts[ftype] += 1
			if fname is not None:
				name_type[fname] = ftype
			if "Parent" in feature.attr:
				features_with_parents.append(feature)
			else:
				feature_types[ftype]["parent unspecified"] += 1
	
	printer.write("Sorting parents...")
	c = 0
	for feature in features_with_parents:
		if c % 10000 == 0:
			printer.write("Processed %s parents..." % c)
		c += 1
		pnames = feature.attr["Parent"]
		ftype = feature.attr["type"]
		if pnames == "":
			feature_types[ftype]["parent unspecified"] += 1
		else:
			if len(pnames) > 1:
				feature_types[ftype]["multiple parents"] += 1
			else:
				ptype = name_type.get(pnames[0],"parent not in database")
				feature_types[ftype][ptype] += 1

	rows = sorted(feature_types.keys())
	cols = rows + ["parent unspecified","parent not in database","multiple parents"]

	with argsopener(args.outfile,args,"w") as fh:
		printer.write("Writing %s..." % args.outfile)
		header = "#feature_type\tcount\t" + "\t".join(cols) + "\n"
		fh.write(header)
		for r in rows:
			sout = "%s\t%s" % (r, feature_counts[r])
			for i in cols:
				sout += "\t%s" % feature_types[r].get(i,0)
			fh.write("%s\n" % sout)

	printer.write("Done.")
	
if __name__ == "__main__":
		main()
