# yeti v0.8.2
(c) 2013 Joshua G Dunn


## Introduction

`yeti` is intended to be:


    1. A lightweight library that makes it easy for biologists and 
       computational biologists to manipulate genomic information using
       familiar tools an idioms.
       
       To this end, `yeti` defines simple APIs that serve as a
       bridge between the rich scientific environment of the SciPy stack,
       and a fair chunk of the variety of file formats that encode genomics
       data (e.g. BAM, Wiggle, bedGraph, GTF2, GFF3, BED, BigBed, PSL, et c).

    2. A set of scripts, built upon this library, that implement common
       sequencing workflows end-to-end. Biologists can use these scripts
       directly on their aligned sequencing data, without much preprocessing
       and certainly without neeing to delve into the APIs described above.



## INSTALLATION

### This isn't presently true
Stable versions can be fetched from [PyPi](https://pypi.python.org)
using `pip` or `easy_install`. Development versions can be fetched
from our [github page]().


## REQUIREMENTS

The following packages and their dependencies are required:

- Python     >= 2.7 or >= 3.3
- Numpy      >= 1.9.2
- Scipy      >= 0.15.1
- Matplotlib >= 1.3
- Pandas     >= 0.16.0
- Pysam      >= 0.7.7
- BioPython  >= 1.64


## DOCUMENTATION

At this point, most of the documentation is in the docstrings.

For info on how to use library components, see modules docstrings and additionally
the command-line scripts in `yeti.bin`

Otherwise, the package is structured as follows:

Package            | Contains
-------------------|--------------------------------------------------------------
yeti        | 
	bin			   | Command-line scripts
	genomics       | Library components specific to genomics
		readers	   | Parsers/readers/writers of various genomics formats (e.g. BED, GTF)
		scriptlib  | Library components for command-line scripts
	util           | Library components 
		io         | Parsers/readers for io streams, writers for files 
		services   | Various other low-level stuff
	test           | Unit tests, integrative tests, & demo usage



## CONTRIBUTING

We welcome your contributions. Please see the documentation for details
on how you can join in.

## LICENSE

This software is released under the BSD 3-Clause license.

## DISCLAIMER

This software is free, not guaranteed, et c et c


1. We follow Vincent Driessen's git branching model, described
[here](http://nvie.com/posts/a-successful-git-branching-model/).
[git-flow](https://github.com/nvie/gitflow/wiki/Installation)
provides tools useful for this.
       
2. Code should be formatted as described
in [PEP8](https://www.python.org/dev/peps/pep-0008),
especially noting that we use four spaces for indentation.

3. Docstrings should be formatted as described using
the [numpydoc plugin](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
for [Sphinx](http://sphinx-doc.org/),
with guidance from [PEP257](https://www.python.org/dev/peps/pep-0257).
This means that docstrings are formatted in
[reStructuredText](http://docutils.sourceforge.net/docs/user/rst/quickref.html).
Beyond that, please reference classes in docstrings using the shortcut `|ClassName|`
rather than `:py:class:\`package.module.ClassName\``, because we perform these substitutions
automatically and generate cross-references in the doc build process. It also makes
docstrings more readable.
