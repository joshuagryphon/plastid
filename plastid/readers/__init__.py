#!/usr/bin/env python
"""
Package overview
================

This package contains parsers for the following file types:

    ======================================    ===================================
    **Annotation files**
    -----------------------------------------------------------------------------
    :py:mod:`plastid.readers.bed`             `BED`_ and `BED X+Y <BED extended format>`_ formats
    :py:mod:`plastid.readers.bigbed`          `BigBed`_                      
    :py:mod:`plastid.readers.gff`             `GTF2`_ and `GFF3`_            
    :py:mod:`plastid.readers.psl`             `PSL`_ (BLAT format)           
    --------------------------------------    -----------------------------------
    **Alignment & count files**     
    -----------------------------------------------------------------------------
    :py:mod:`plastid.readers.bowtie`          `bowtie`_ alignment            
    :py:mod:`plastid.readers.wiggle`          `bedGraph`_ and `Wiggle`_      
    --------------------------------------    -----------------------------------
    **Miscellaneous**
    -----------------------------------------------------------------------------
    :py:mod:`plastid.readers.autosql`          `autoSql`_                     
    ======================================    ===================================


 .. note::

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


Helper code can be found in the following modules:

    =================================    ==========================================
    **Module**                           **Contents**
    ---------------------------------    ------------------------------------------
    :mod:`plastid.readers.common`        Helper functions used by many annotation 
                                         file readers
    
    :mod:`plastid.readers.gff_tokens`    Functions for parsing and exporting
                                         attributes in the ninth column of 
                                         `GTF2`_ and `GFF3`_ files
    =================================    ==========================================
"""
