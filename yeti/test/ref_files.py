#!/usr/bin/env python
"""This module contains paths to reference files used by various unit and funcitonal tests
"""
import os
from pkg_resources import resource_filename

#===============================================================================
# INDEX: small data set for faster/unit tests
#===============================================================================

MINI_PATH = resource_filename("plastid","test/data/mini")
MINI = {
    'vec_center_0_rc.txt'  : 'count_vectors/center_0_rc.txt',
    'vec_center_12_rc.txt' : 'count_vectors/center_12_rc.txt',
    'vec_center_0_fw.txt'  : 'count_vectors/center_0_fw.txt',
    'vec_center_12_fw.txt' : 'count_vectors/center_12_fw.txt',

    'vec_fiveprime_0_fw.txt'  : 'count_vectors/fiveprime_0_fw.txt',
    'vec_fiveprime_0_rc.txt'  : 'count_vectors/fiveprime_0_rc.txt',
    'vec_fiveprime_15_fw.txt' : 'count_vectors/fiveprime_15_fw.txt',
    'vec_fiveprime_15_rc.txt' : 'count_vectors/fiveprime_15_rc.txt',

    'vec_threeprime_0_rc.txt'  : 'count_vectors/threeprime_0_rc.txt',
    'vec_threeprime_0_fw.txt'  : 'count_vectors/threeprime_0_fw.txt',
    'vec_threeprime_15_fw.txt' : 'count_vectors/threeprime_15_fw.txt',
    'vec_threeprime_15_rc.txt' : 'count_vectors/threeprime_15_rc.txt',

    'fasta_file'          : 'fasta/chrA.fa',
    'bamfile'             : 'align/chrA_tophat.bam',
    'bowtie_file'         : 'align/chrA_unspliced.bowtie',
    'bed_file'            : 'bed/chrA.bed',
    'bigbed_file'         : 'bed/chrA.bb',
    'gtf2_file'           : 'bed/chrA.gtf',
    'gff3_file'           : 'bed/chrA.gff',
    'psl_file'            : 'bed/chrA.psl',

    'wig_bedgraph'         : 'wig/bedgraph',
    'wig_bedgraph_fw'      : 'wig/bedgraph_fw.wig',
    'wig_bedgraph_rc'      : 'wig/bedgraph_rc.wig',
    'wig_variable_step'    : 'wig/variable_step',
    'wig_variable_step_fw' : 'wig/variable_step_fw.wig',
    'wig_variable_step_rc' : 'wig/variable_step_rc.wig',
}

MINI = { K : os.path.join(MINI_PATH,V) for K,V in MINI.items() }
"""Dictionary of paths to files describing a small set of 28-mer alignments"""


#===============================================================================
# INDEX: big data set for functional tests
#===============================================================================

RPATH = resource_filename("plastid","test/data")

REF_FILES = { 
             # full yeast annotation
              "yeast_gtf"           : os.path.join(RPATH,"annotations","sgd_plus_utrs.gtf"),
              "yeast_gff"           : os.path.join(RPATH,"annotations","sgd_plus_utrs.gff"),

              # 100 yeast transcripts, used in testing various readers
              "100transcripts_gtf" :  os.path.join(RPATH,"annotations","100transcripts.gtf"),
              "100transcripts_gff" :  os.path.join(RPATH,"annotations","100transcripts.gff"),
              "100transcripts_bed" :  os.path.join(RPATH,"annotations","100transcripts.bed"),
              "100cds_bed"          :  os.path.join(RPATH,"annotations","100cds.bed"),
              "100cds_antisense_bed":  os.path.join(RPATH,"annotations","100cds_antisense.bed"),

              "100transcripts_bigbed" :  os.path.join(RPATH,"annotations","100transcripts.bb"),
              "100cds_bigbed"          :  os.path.join(RPATH,"annotations","100cds.bb"),
              "100cds_antisense_bigbed":  os.path.join(RPATH,"annotations","100cds_antisense.bb"),
              
              "100transcripts_gtf_tabix" : os.path.join(RPATH,"annotations","100transcripts_sorted.gtf.gz"),
              
              "100transcripts_bed_tabix" : os.path.join(RPATH,"annotations","100transcripts_sorted.bed.gz"),
              "100cds_bed_tabix" : os.path.join(RPATH,"annotations","100cds_sorted.bed.gz"),
              "100as_cds_bed_tabix" : os.path.join(RPATH,"annotations","100cds_antisense_sorted.bed.gz"),

              # reference output files from command-line scripts
              "yeast_parent_child"  : os.path.join(RPATH,"annotations","sgd_plus_utrs_gff3_parent_child_table.txt"),
              "yeast_parent_child_exclude_cols"  : os.path.join(RPATH,"annotations","sgd_plus_utrs_gff3_parent_child_table_exclude.txt"),
              "yeast_fasta"         : os.path.join(RPATH,"annotations","2013-01-23.saccharomyces_cerevisiae.fa"),
              "yeast_twobit"        : os.path.join(RPATH,"annotations","2013-01-23.saccharomyces_cerevisiae.2bit"),

              "yeast_juncs"         : os.path.join(RPATH,"command_line","sgd_plus_utrs.juncs"),
              "yeast_crossmap_o12_26_0"    : os.path.join(RPATH,"command_line","2013-01-23.sc_o12_26_0_crossmap.bed"),
              "yeast_crossmap_o12_26_0_bb" : os.path.join(RPATH,"command_line","2013-01-23.sc_o12_26_0_crossmap_sorted.bb"),

              "yeast_crossmap_o0_26_0"    : os.path.join(RPATH,"command_line","2013-01-23.sc_o0_26_0_crossmap.bed"),
              "yeast_crossmap_o12_26_2"   : os.path.join(RPATH,"command_line","2013-01-23.sc_o12_26_2_crossmap.bed"),
              
              "yeast_rp_bam"        : os.path.join(RPATH,"command_line","gen_reads.bam"),
              "yeast_psite"         : os.path.join(RPATH,"command_line","gen_p_offsets.txt"),
              "yeast_mini_bed"      : os.path.join(RPATH,"command_line","transcripts_chosen.bed"),
              "yeast_phasing"       : os.path.join(RPATH,"command_line","gen_phasing.txt"),
              
              "yeast_metagene_cds_start" : os.path.join(RPATH,"command_line","gen_cds_start_rois.txt"),
              "yeast_metagene_cds_stop"  : os.path.join(RPATH,"command_line","gen_cds_stop_rois.txt"),              
              "yeast_metagene_cds_start_bed" : os.path.join(RPATH,"command_line","gen_cds_start_rois.bed"),
              "yeast_metagene_cds_stop_bed"  : os.path.join(RPATH,"command_line","gen_cds_stop_rois.bed"),
              
              "yeast_metagene_cds_start_profile" : os.path.join(RPATH,"command_line","gen_cds_start_metagene_profile.txt"),
              "yeast_metagene_cds_stop_profile"  : os.path.join(RPATH,"command_line","gen_cds_stop_metagene_profile.txt"),
              
              "yeast_metagene_cds_start_normcounts" : os.path.join(RPATH,"command_line","gen_cds_start_normcounts.txt"),
              "yeast_metagene_cds_stop_normcounts"  : os.path.join(RPATH,"command_line","gen_cds_stop_normcounts.txt"),
              
              "yeast_metagene_cds_start_rawcounts" : os.path.join(RPATH,"command_line","gen_cds_start_rawcounts.txt"),
              "yeast_metagene_cds_stop_rawcounts"  : os.path.join(RPATH,"command_line","gen_cds_stop_rawcounts.txt"),

              "yeast_findjuncs_bed"   :  os.path.join(RPATH,"annotations","sgd_plus_utrs_findjuncs.bed"),
              "yeast_findjuncs_juncs" :  os.path.join(RPATH,"annotations","sgd_plus_utrs_findjuncs.juncs"),
              "yeast_findjuncs_top"   :  os.path.join(RPATH,"annotations","sgd_plus_utrs_findjuncs_top250.bed"),
              "yeast_findjuncs_bot" :  os.path.join(RPATH,"annotations","sgd_plus_utrs_findjuncs_bot230.bed"),

              "yeast_counts_in_region_mask" : os.path.join(RPATH,"command_line","gen_counts_in_region_mask.txt"),
              "yeast_counts_in_region_no_mask" : os.path.join(RPATH,"command_line","gen_counts_in_region_no_mask.txt"),

              'slidejuncs_seqs'                        : os.path.join(RPATH,'annotations','juncs','junction_seqs.fa'),
              'slidejuncs_input'                       : os.path.join(RPATH,'annotations','juncs','input.bed'),
              'slidejuncs_ref'                         : os.path.join(RPATH,'annotations','juncs','all_known_ref.bed'),
              'slidejuncs_crossmap'                    : os.path.join(RPATH,'annotations','juncs','junc_crossmap.bed'),
              'slidejuncs_known_juncs_non_crossmap'    : os.path.join(RPATH,'annotations','juncs','known_juncs_non_crossmapping.bed'),
              'slidejuncs_known_juncs_crossmap'        : os.path.join(RPATH,'annotations','juncs','known_juncs_crossmapping.bed'),
              'slidejuncs_to_slide_known_non_crossmap' : os.path.join(RPATH,'annotations','juncs','to_slide_non_crossmapping.bed'),
              'slidejuncs_to_slide_known_crossmap'     : os.path.join(RPATH,'annotations','juncs','to_slide_crossmapping.bed'),
              'slidejuncs_noncan_no_ref'               : os.path.join(RPATH,'annotations','juncs','noncan_no_ref.bed'),
              'slidejuncs_expected_untouched'          : os.path.join(RPATH,'annotations','juncs','expected_untouched.bed'),
    
            }

COUNT_OPTIONS = " --count_files %s --fiveprime_variable --offset %s" % (REF_FILES["yeast_rp_bam"],REF_FILES["yeast_psite"])
ANNOTATION_OPTIONS = " --annotation_format GTF2 --add_three --annotation_files %s" % REF_FILES["yeast_gtf"]
MASK_OPTIONS = " --mask_annotation_format BED --mask_annotation_files %s" % REF_FILES["yeast_crossmap_o12_26_0"]
