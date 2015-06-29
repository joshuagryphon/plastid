#!/usr/bin/env python
"""
Package overview
================

This package contains parsers for various file types, such as:

    =========================    ==============================================
    **Annotation files**
    ---------------------------------------------------------------------------
    `BED`_ and BED+x formats     :py:mod:`yeti.readers.bed`
    `BigBed`_                    :py:mod:`yeti.readers.bigbed`
    `GTF2`_ and `GFF3`_          :py:mod:`yeti.readers.gff`
    `PSL`_ (BLAT format)         :py:mod:`yeti.readers.psl`
    -------------------------    ----------------------------------------------
    **Alignment & count files**     
    ---------------------------------------------------------------------------
    `bowtie`_ alignment          :py:mod:`yeti.readers.bowtie`
    `bedGraph`_ and `Wiggle`_    :py:mod:`yeti.readers.wiggle`
    -------------------------    ----------------------------------------------
    **Miscellaneous**
    ---------------------------------------------------------------------------
    `autoSql`_                   :py:mod:`yeti.readers.autosql`
    =========================    ==============================================


In addition:

  - `BAM`_ files are supported via `Pysam`_ and may be accessed through |BAMGenomeArray|.
  - All of the readers above (except autoSql) are compatible with `tabix`_
    compression, also supported via `Pysam`_. To use one of these with
    a `tabix`_-compressed file, wrap the open filehandle with
    :py:obj:`pysam.tabix_iterator`, and pass a dummy function to
    retrieve the uncompressed text::

        dummy_func = lambda x, y : x
        my_reader = BED_Reader(pysam.tabix_iterator(open("some_file.bed",dummy_func)))


"""
