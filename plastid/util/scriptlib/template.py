#!/usr/bin/env python
"""This is a template for a command-line script. To write a script, edit this file,
making changes as necessary, Make changes below as necessary, and, when finished,
update this docstring.
"""
import inspect
import sys
import argparse
from plastid.util.scriptlib.help_formatters import format_module_docstring
from plastid.util.io.filters import NameDateWriter
from plastid.util.io.openers import get_short_name

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


def main(argv=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
    argv : list, optional
        A list of command-line arguments, which will be processed
        as if the script were called from the command line if
        :py:func:`main` is called directly.

        Default: sys.argv[1:] (actual command-line arguments)
    """
    # fill out command-line program here
    
    # add parents from plastid.scriptlib.argparsers as necessary
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[])

    # add your own argument
    parser.add_argument("--foo",type=str,help="Some argument foo")

    args = parser.parse_args(argv)
    
    # write program body
    pass


if __name__ == "__main__":
    main()