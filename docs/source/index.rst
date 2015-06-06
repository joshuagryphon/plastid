.. yeti documentation master file, created by
   sphinx-quickstart on Fri Dec  5 11:55:54 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

yeti v0.5 welcome!
==========================

Introduction
------------

:py:data:`yeti` is a lightweight library designed to facilitate analysis
of high throughput sequencing data, especially, but not limited to, ribosome profiling.

The primary design goal is to convert genomic data into Pythonic objects that can be
manipulated by preexisting tools in the `SciPy stack <http://www.scipy.org/stackspec.html>`_,
introducing the smallest possible number of new classes or data types.

To this end, :py:data:`yeti` provides:

	* :doc:`command-line scripts </cli_howto>` that implement common sequencing workflows,
	  as well as several analyses specific to ribosome profiling
	
	* A library of :ref:`data structures and methods <overview-of-data-structures>`
	  for interactive or ad-hoc analyses, 
	  as well as :doc:`readers for various file formats </generated/yeti.readers>`
	  (`Wiggle`_, `bedGraph`_, `bowtie`_, `BED`_, `BigBed`_, `GTF2`_, `GFF3`_, `PSL`_; `BAM`_ supported via `Pysam`_).
	 
	* Components to simplify :doc:`writing your own command-line scripts </generated/yeti.util.scriptlib>` 
	  
	  
For more info, see:

	* :doc:`Getting started </tour>` for an overview of some of the most useful data structures
	
	* :doc:`Walk-throughs </interactive_howto>` that guide you through several interactive analyses
	
	* Brief descriptions of the :doc:`command-line scripts </cli_howto>` that are included

	* Complete :ref:`API documentation <modindex>` for developers

   
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

