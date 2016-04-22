#!/usr/bin/env python
"""
Package overview
================

This package contains parsers for various file types. All parsers behave as
iterators, and return data
that is standarized to 0-indexed, half-open coordinate systems, in keeping with
Python  conventions. In addition, all of the annotation file readers
(except `BigBed`_) are compatible with `tabix`_ compression, which is supported
via `Pysam`_. 


    ======================================    =======================================
    **Annotation files**
    ---------------------------------------------------------------------------------
    :py:mod:`plastid.readers.bed`             `BED`_ and :term:`Extended BED` formats
    :py:mod:`plastid.readers.bigbed`          `BigBed`_                      
    :py:mod:`plastid.readers.gff`             `GTF2`_ and `GFF3`_            
    :py:mod:`plastid.readers.psl`             `PSL`_ (BLAT format)           
    --------------------------------------    ---------------------------------------
    **Alignment & count files**     
    ---------------------------------------------------------------------------------
    :py:mod:`plastid.readers.bowtie`          `bowtie`_ alignment            
    :py:mod:`plastid.readers.wiggle`          `bedGraph`_ and `Wiggle`_
    :py:mod:`plastid.readers.bigwig`          `BigWig`_
    --------------------------------------    ---------------------------------------
    **Miscellaneous**
    ---------------------------------------------------------------------------------
    :py:mod:`plastid.readers.autosql`          `autoSql`_                     
    ======================================    =======================================



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
