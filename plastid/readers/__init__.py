#!/usr/bin/env python
"""
Package overview
================

This package contains parsers for the following file types:

    =========================    ==============================================
    **Annotation files**
    ---------------------------------------------------------------------------
    `BED`_ 4-12 formats          :py:mod:`plastid.readers.bed`
    `BigBed`_                    :py:mod:`plastid.readers.bigbed`
    `GTF2`_ and `GFF3`_          :py:mod:`plastid.readers.gff`
    `PSL`_ (BLAT format)         :py:mod:`plastid.readers.psl`
    -------------------------    ----------------------------------------------
    **Alignment & count files**     
    ---------------------------------------------------------------------------
    `bowtie`_ alignment          :py:mod:`plastid.readers.bowtie`
    `bedGraph`_ and `Wiggle`_    :py:mod:`plastid.readers.wiggle`
    -------------------------    ----------------------------------------------
    **Miscellaneous**
    ---------------------------------------------------------------------------
    `autoSql`_                   :py:mod:`plastid.readers.autosql`
    =========================    ==============================================


Helper code can be found in the following modules:

    ==============================    ==========================================
    **Module**                        **Contents**
    ------------------------------    ------------------------------------------
    :mod:`plastid.readers.common`        Helper functions used by many annotation 
                                      file readers
    
    :mod:`plastid.readers.gff_tokens`    Functions for parsing and exporting
                                      attributes in the ninth column of 
                                      `GTF2`_ and `GFF3`_ files
    ==============================    ==========================================

Notes
-----
  
  - All readers return objects whose genome coordinates are *0-indexed and
    half-open,* in keeping with Python conventions, *regardless* of the representation
    in the underlying file.
    
  - All of the annotation file readers (except `BigBed`_) are compatible with `tabix`_
    compression (supported via `Pysam`_). To use one of these with
    a `tabix`_-compressed file, pass the keyword argument `tabix=True`,
    wrap the open filehandle with :py:obj:`pysam.tabix_iterator`::

        >>> import pysam
        >>> my_reader = BED_Reader(pysam.tabix_iterator(open("some_file.bed"),pysam.asTuple()),tabix=True)

  - `BAM`_ files are supported via `Pysam`_ and may be accessed through |BAMGenomeArray|.

"""
