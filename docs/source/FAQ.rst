Frequently asked questions
==========================

.. contents:: FAQ contents
    :depth: 3
    :local:


.. _faq-run:

Installation and runtime
------------------------


.. _faq-install-fails:

Installation fails in pip with no obvious error message
.......................................................

Installation can fail for multiple reasons. To figure out what is responsible,
repeat installation passing the ``--verbose`` flag to ``pip``: 

.. code-block:: shell

   $ pip install --no-cache-dir --verbose plastid | tee 2>&1 plastid_install_log.txt

Then find the corresponding error message below. If the error is not listed,
let us know by filing a bug report at `our issue tracker`_. Please attach
`plastid_install_log.txt` to your report to help us figure out what is going
on.


.. _faq-install-numpy-first:

Installer quits with an error message about `numpy`_, `pysam`_ and `Cython`_
............................................................................

The most common cause of installation failure is that the setup script for
`plastid` requires `numpy`_, `pysam`_, and `Cython`_ to be installed before it
can be run. In this case, the following message should appear:

.. code-block:: none

   *** IMPORTANT INSTALLATION INFORMATION ***

   plastid setup requires numpy>=1.9.0, pysam>=0.8.4, and cython>=0.22 to be preinstalled. Please
   install these via pip, and retry:

       $ pip install --upgrade numpy pysam cython
       $ pip install plastid

If this is the problem, simply install `numpy`_, `pysam`_ and `Cython`_ first,
then repeat the installation:

.. code-block:: shell

   $ pip install numpy pysam cython
   $ pip install plastid

Check if install worked:

.. code-block:: shell
 
   $ pip list
    
You should see `numpy`_, `pysam`_, and `Cython`_ in the list.


.. _faq-install-fails-with-prereqs:

`numpy`_, `pysam`_ and `Cython`_ are up to date, but install still fails
........................................................................

If install fails when `numpy`_, `pysam`_ and `Cython`_ are already up to date,
then there may be multiple versions of these libraries installed on your system,
with :data:`plastid` seeing an earlier version than `pip`. A few users have
come across this problem when installing :data:`plastid` in a `conda`_/`Anaconda`_
environment.

A solution is to install inside a vanilla environment in a fresh `virtualenv`_
(note: *not* a `conda`_/`Anaconda`_ environment). See instructions at
:ref:`Install inside a virtualenv <install-inside-venv>`.

If install succeeds in a `virtualenv`_, this suggests that there are in fact
multiple versions of one or more of plastid's dependencies installed on your
system. In this case, ``plastid`` can be used inside the `virtualenv`_.



.. _faq-install-conda:

Install fails inside a `conda`_/`Anaconda`_ environment
.......................................................

Two users have reported difficulties installing inside `conda`_/`Anaconda`_
environments. Despite having up-to-date versions of `numpy`_, `pysam`_, and
`Cython`_ installed in the `conda`_ environment, during the build process
``pip`` found an incompatible version of `Cython`_ or of `pysam`_.

A workaround is to install inside a `virtualenv`_, instructions for which can be
found at :ref:`Install inside a virtualenv <install-inside-venv>`.


.. _faq-locale-error-osx:
 
Locale error when installing or running ``plastid`` on OSX
..........................................................

This is known to occur on OSX. In this case, you should see a stack trace ending
with something like:

.. code-block:: none
 
   from docutils.utils.error_reporting import locale_encoding, ErrorString, ErrorOutput
     File "/Applications/anaconda/lib/python2.7/site-packages/docutils/utils/error_reporting.py", line 47, in <module>
       locale_encoding = locale.getlocale()[1] or locale.getdefaultlocale()[1]
     File "/Applications/anaconda/lib/python2.7/locale.py", line 543, in getdefaultlocale
       return _parse_localename(localename)
     File "/Applications/anaconda/lib/python2.7/locale.py", line 475, in _parse_localename
       raise ValueError, 'unknown locale: %s' % localename
   ValueError: unknown locale: UTF-8

Please see the workaround 
`here <http://blog.remibergsma.com/2012/07/10/setting-locales-correctly-on-mac-osx-terminal-application/>`_.


.. _faq-macintosh-cflags:
 
Install fails on OSX with `error code 1`
........................................

If installing on OSX and you find an error message that resembles the following:

.. code-block:: none
 
   Command "/usr/local/opt/python/bin/python2.7 -c "import setuptools, tokenize;\
   __file__='/private/var/folders/8y/xm0qbq655f1d4v20kq5vvfgm0000gq/T/pip-build-0bVdPy/pysam/setup.py';\
   exec(compile(getattr(tokenize, 'open', open)(__file__).read().replace('\r\n', '\n'), __file__, 'exec'))"\
   
    install --record /var/folders/some-folder/install-record.txt --single-version-externally-managed \
            --compile --user --prefix=" failed with error code 1 in /private/var/folders/some-folder/pysam

Before installing, type:

.. code-block:: shell
 
   $ export CFLAGS=-Qunused-arguments
   $ export CPPFLAGS=-Qunused-arguments

and then retry.


.. _faq-script-path:
 
``command not found``: I can't run any of the command line scripts
..................................................................

If you receive a `command not found` error from the shell, the folder containing
the command-line scripts might not be in your environment's ``PATH`` variable.

Command-line scripts will be installed wherever your system configuration dictates.
On OSX and many varities of linux, the install path for a single-user install is
``~/bin`` or ``~/.local/bin``. For system-wide installs, the path is typically
``/usr/local/bin``. Make sure the appropriate location is in your ``PATH`` by
the following line adding to your ``.bashrc``, ``.bash_profile``, or ``.profile``
(depending on which your system uses):

.. code-block:: shell

   export PATH=~/bin:~/.local.bin:/usr/local/bin:$PATH


.. _faq-script-insufficient-arguments:
 
A script won't run, reporting `error: too few arguments`
........................................................

If you see the following error:

.. code-block:: none
 
   <script name>: error: too few arguments

Try re-ordering the script arguments, so that all of the required arguments
(the ones that don't start with ``-``) come first. For example, change:

.. code-block:: shell
 
   $ cs count --fiveprime --offset 13 --min_length 23 --max_length 35 \
              --count_files ../some_file.bam some_file.positions some_sample_name
 
to

.. code-block:: shell

   $ cs count some_file.positions some_sample_name \
              --fiveprime --offset 13 --min_length 23 --max_length 35 \
              --count_files ../some_file.bam 

Alternatively, put a ``--`` before the required options:

.. code-block:: shell

   $ cs count --fiveprime --offset 13 --min_length 23 --max_length 35 \
              --count_files ../some_file.bam \
              -- some_file.positions some_sample_name

Things should then run.



.. _faq-distribution-error: 

I get an ``ImportError`` or ``DistributionError`` when using :data:`plastid`
............................................................................

If you get an error like the following::

    Traceback (most recent call last):
       File "/home/user/Rib_prof/venv/bin/crossmap", line 5, in <module>
         from pkg_resources import load_entry_point
       File "/home/user/Rib_prof/venv/lib/python2.7/site-packages/pkg_resources/__init__.py", line 2970, in <module>
         working_set = WorkingSet._build_master()
       File "/home/user/Rib_prof/venv/lib/python2.7/site-packages/pkg_resources/__init__.py", line 567, in _build_master
         ws.require(__requires__)
       File "/home/user/Rib_prof/venv/lib/python2.7/site-packages/pkg_resources/__init__.py", line 876, in require
         needed = self.resolve(parse_requirements(requirements))
       File "/home/user/Rib_prof/venv/lib/python2.7/site-packages/pkg_resources/__init__.py", line 761, in resolve
         raise DistributionNotFound(req)
     pkg_resources.DistributionNotFound: scipy>=0.12.0 


One or more dependencies (in this example, `SciPy`_ is not installed).
Please see :ref:`faq-install-numpy-first`, above.




-------------------------------------------------------------------------------

.. _faq-analysis:
 
Analysis
--------

.. _faq-dutp-data:
 
Can I use ``plastid`` with reverse-complemented sequencing data, like dUTP sequencing?
......................................................................................

Yes.

Kits like Illumina's `Truseq Stranded mRNA Library Prep Kit <http://www.illumina.com/products/truseq_stranded_mrna_library_prep_kit.html>`_,
yield reads that are anti-sense to the mRNA from which they were generated, so
the data coming off the sequencer will be reverse-complemented compared to the
original strand that was cloned.

To use this data in :data:`plastid`, reverse-complement your FASTQ file using the
`fastx_reverse_complement <http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_revcomp_usage>`_
tool from the Hanon lab's `fastx toolkit`_. Then align the reverse-complemented
data using your favorite aligner.


.. _faq-paired-end:

Can ``plastid`` be used with paired-end data?
.............................................

Yes, but two points:

 - Because there are few nucleotide-resolution assays that used paried-end sequencing,
   it has been unclear what sorts of :term:`mapping functions <mapping funcgtion>`
   might be useful.
   
   If you have a suggestion for one, please submit your suggestion
   with a use case on `our issue tracker`_. It would be helpful!

 - For simple gene expression counting, it is possible to use the ``--fiveprime`` 
   (implemented in :class:`~plastid.genomics.roitools.FivePrimeMapFactory`)
   mapping function with zero offset. Accuracy can be improved by counting a single
   read from each pair in which both reads are mapped.
   
   However, it is critical to retain the read that appears on the same strand
   as the gene from which it arose.
   
   If the library was prepared using dUTP chemistry, as in many paired-end 
   prep kits, select `read2` from each pair using ``samtools``:
   
   .. code-block:: shell

      $ samtools view  -f 129 -b -o single_from_pair.bam paired_end_file.bam
       
   Otherwise, select `read1`:
    
   .. code-block:: shell

      $ samtools view  -f 65 -b -o single_from_pair.bam paired_end_file.bam
     
  Then use `single_from_pair.bam` with ``plastid`` as usual.



.. _faq-psite-use-aggregate:

The P-site and/or Metagene scripts show few or zero read in their output
........................................................................

This occurs in datasets with few counts, because |psite| and |metagene| plots
the median density at each position. In this case, there are a few options:

 - increase the minimum counts required to be included in the metagene / P-site
   estimate. Set ``--min_counts`` argument to a high number (e.g. for a 100 nt
   normalization region, choose >= 100 counts)
 
 - the metagene profile or P-site can be estimated from aggregate counts (as
   opposed to median density) at each position using the ``--aggregate`` argument,
   as shown :ref:`here <psite-use-aggregate>`.
   
   This might add some noise to the data, but it should still be interpretable
   if the gene models are good



.. _faq-zero-counts-because-antisense:
 
``cs``, ``counts_in_region``, or some other part of ``plastid`` reports zero counts for my gene, even though there are read alignments there
............................................................................................................................................

The default behavior for all of the scripts and tools in :data:`plastid` is to 
exclude reads that are antisense to any given genomic feature when calculating
coverage over that feature.

Paired end libraries, and single-end libraries that have been prepared with dUTP
sequencing or a number of other protocols will contain read alignments antisense
to the original mRNAs, causing these reads to be considered antisense to genes,
and therefore excluded from gene expression totals.

See :ref:`faq-dutp-data` for instructions on how to reverse-complement your data
for single-end dUTP data; or :ref:`faq-paired-end` for info on using 
:data:`plastid` with paired-end data.


.. _faq-analysis-fractional-counts:

Why do some scripts report fractional count numbers?
....................................................

Fractional counts for :term:`read alignments` arise when using a alignment
:term:`mapping rule` that maps reads fractionally over multiple positions (such
as ``--center`` mapping). See the  discussion of :doc:`concepts/mapping_rules`,
where these are discussed in depth.




.. _faq-igv-vs-mapped-wiggle:

Why does `IGV`_ report higher coverage at a given nucleotide than the file exported from ``make_wiggle``?
.........................................................................................................
When `IGV`_ calculates coverage of a nucleotide, it counts the number of alignments covering that nucleotide.
So, a 30-nucleotide read would contribute 30 :term:`counts` to a dataset.

While it is possible to write any mapping rule in :mod:`plastid`, the :term:`mapping rules <mapping rule>`
included by default count each read only once (e.g. at their 5' end, 3' end, et c). Even when using
*center* or *entire* mapping, each position covered by a read alignment is only incremented by :math:`1.0/\ell`,
where :math:`\ell` is the length of the read. So, in this case, a 30-nucleotide read would only 
contribute 1 :term:`count <counts>` to a dataset. See :doc:`/concepts/mapping_rules/` for more information.


.. _faq-cs-vs-counts-in-region:

What are the differences between ``counts_in_region`` and ``cs``?
.................................................................
:mod:`~plastid.bin.counts_in_region` very simply counts read coverage (or any data) over
regions of interest, and reports those numbers in terms of :term:`counts` and :term:`RPKM`. It can 
optionally take a :term:`mask file`, if there are genomic positions in the regions
of interest which should be excluded from analysis. Otherwise, it makes no corrections.

:mod:`~plastid.bin.cs` is more complex, and is principally designed to make rough estimates
of gene expression at the gene, rather than transcript, level. In so doing, it makes several
heuristic corrections to regions before tabulating their :term:`counts` and :term:`RPKM`. Specifically:

 #. Genes that have transcripts that share exons are merged into single entities

 #. Gene areas are defined for each merged geen by including all positions occupied
    by all transcripts from that merged gene

 #. Regions occupied by two or more merged genes on the same strand are excluded from
    the calculation of expression values for both genes
 
 #. Optionally, a :term:`mask file` can be used to exclude any other positions from
    analysis.

 #. Expression values (in :term:`counts` and :term:`RPKM`) are tabulated for the entire
    gene area (reported as *exon_counts* and *exon_rpkm*) as well as for sub regions,
    if the gene is coding. Specifically, *cds_counts* and *cds_rpkm* are calculated
    from counts that cover positions in the gene area that are annotated as CDS in
    **all** transcripts in the merged gene. Ditto for 5' and 3' UTRs

Either one can be an appropriate starting place for a pipeline, depending upon your needs.
See the documentation and/or source code for |cs| and |counts_in_region| for further
discussion. 


.. _faq-segmentchain-gff3:

Why does ``SegmentChain.as_gff3()`` sometimes throw errors?
...........................................................

The incredible fle flexibility of the `GFF3`_ file format introduces ambiguities
for representation of discontinuous features: some sort of parent-child relationship
needs to exist, and, except in the case of transcripts, :data:`plastid` doesn't
know which one to use.  

See :ref:`this advice <data-export-gff3>` on how to handle this.


.. _faq-analysis-deseq:

How do I prepare data for differential gene expression analysis in `DESeq`_?
............................................................................

See :doc:`examples/gene_expression` in the :doc:`examples` section.




-------------------------------------------------------------------------------

.. _faq-tests:

Tests
-----

The tests won't run
...................
In order to run the tests, you need to download the `test dataset`_ and unpack
it into ``plastid/test/``.

We decided not to include the test data in the main package in order to keep
the package download and the github repository small.


-------------------------------------------------------------------------------


.. toctree::
   :maxdepth: 2
