TABLE?=gene_count_table.txt
GTF=../annotations/sgd_plus_utrs.gtf
CROSSMAP=2013-01-23.sc_o12_26_0_crossmap.bed
CROSSMAP_FORMAT=BED

# Phony targets
all : gen_p_offsets.txt \
	  gen_cds_start_metagene_profile.txt \
	  gen_cds_stop_metagene_profile.txt \
	  gen_cds_start_metagene.svg \
	  gen_cds_stop_metagene.svg \
	  gen_cs_count_gene.txt \
	  gen_reads_fiveprime_variable_fw.wig \
	  gen_reads_fiveprime_variable_rc.wig \
	  gene_count_table.txt \
	  gen_phasing.txt \
      gen_counts_in_region_mask.txt \
      gen_counts_in_region_no_mask.txt

reads : gen_reads.fa

align : gen_reads.bam

counts_in_region : gen_counts_in_region_mask.txt \
                   gen_counts_in_region_no_mask.txt

count_vectors : gen_count_vectors

cs : gen_cs_count_gene.txt

metagene : gen_cds_start_metagene_profile.txt \
		   gen_cds_stop_metagene_profile.txt \
		   gen_cds_start_metagene.svg \
		   gen_cds_stop_metagene.svg 

phasing : gen_phasing.txt

psite : gen_p_offsets.txt

wiggle : gen_reads_fiveprime_variable_fw.wig gen_reads_fiveprime_variable_rc.wig

copy_wiggle :
	cp ../gen*generated*.wig .

clean :
	@rm -rf gen*

.PHONY : clean all reads align count_vectors cs metagene phasing psite wiggle copy_wiggle



# Non-phony targets

gen_reads.fa gene_count_table.txt transcripts_chosen.bed :
	cp ../$@ .

gen_reads.bam : gen_reads.fa
	tophat --bowtie1 --read-mismatches=0 --library-type=fr-firststrand \
		   -p4 --no-novel-juncs \
		   --raw-juncs=$(F250)/annotations/sgd_plus_utrs.juncs \
		   $(F250)/ebwt/2013-01-23.saccharomyces_cerevisiae \
		   $^
	#mv tophat_out/accepted_hits.bam ./gen_reads.bam
	samtools view tophat_out/accepted_hits.bam -H >tmp.sam
	samtools view -F 256 tophat_out/accepted_hits.bam >>tmp.sam
	samtools view -S tmp.sam -b -o $@
	samtools index $@
	#rm -rf tophat_out

gen_cds_start_rois.txt : 
	metagene generate gen_cds_start --downstream 100 \
			 --annotation_format GTF2 \
			 --annotation_files $(GTF) \
			 --add_three \
	         --mask_annotation_files $(CROSSMAP) \
			 --mask_annotation_format $(CROSSMAP_FORMAT)

gen_cds_stop_rois.txt :
	metagene generate gen_cds_stop --upstream 100 \
			 --annotation_format GTF2 \
			 --annotation_files $(GTF) \
			 --add_three \
			 --landmark cds_stop \
	         --mask_annotation_files $(CROSSMAP) \
			 --mask_annotation_format $(CROSSMAP_FORMAT)

gen_p_offsets.txt : gen_reads.bam gen_cds_start_rois.txt
	psite gen_cds_start_rois.txt gen --count_files gen_reads.bam \
		--require_upstream \
		--min_length 26 --max_length 31 \
		--norm_region 70 150
	@rm *normcounts.txt
	@rm *rawcounts.txt

gen_cds_start_metagene_profile.txt : gen_reads.bam gen_p_offsets.txt gen_cds_start_rois.txt
	metagene count gen_cds_start_rois.txt gen_cds_start \
		--count_files gen_reads.bam \
		--fiveprime_variable --offset gen_p_offsets.txt \
		--norm_region 70 150 \

gen_cds_stop_metagene_profile.txt : gen_reads.bam gen_p_offsets.txt gen_cds_stop_rois.txt
	metagene count gen_cds_stop_rois.txt gen_cds_stop \
		--count_files gen_reads.bam \
		--fiveprime_variable --offset gen_p_offsets.txt \
		--norm_region 0 80

gen_cds_start_metagene.svg : gen_cds_start_metagene_profile.txt
	metagene chart $@ $^ --landmark "CDS start" --title "Test metagene at start codon"

gen_cds_stop_metagene.svg : gen_cds_stop_metagene_profile.txt
	metagene chart $@ $^ --landmark "CDS stop" --title "Test metagene at stop codon"

gen_cs_gene.positions :
	cs generate gen_cs --annotation_format GTF2 \
		--annotation_files $(GTF) --add_three \
		--mask_annotation_files $(CROSSMAP) \
		--mask_annotation_format $(CROSSMAP_FORMAT)

gen_cs_count_gene.txt : gen_reads.bam gen_cs_gene.positions gen_p_offsets.txt
	cs count gen_cs_gene.positions gen_cs_count_gene \
		--count_files gen_reads.bam \
		--fiveprime_variable --offset gen_p_offsets.txt \

gen_counts_in_region_no_mask.txt : gen_reads.bam gen_p_offsets.txt transcripts_chosen.bed
	counts_in_region $@ \
		--count_files gen_reads.bam \
		--fiveprime_variable --offset gen_p_offsets.txt \
		--annotation_files transcripts_chosen.bed --annotation_format BED

gen_counts_in_region_mask.txt : gen_reads.bam gen_p_offsets.txt transcripts_chosen.bed
	counts_in_region $@ \
		--count_files gen_reads.bam \
		--fiveprime_variable --offset gen_p_offsets.txt \
		--annotation_files transcripts_chosen.bed --annotation_format BED \
		--mask_annotation_format $(CROSSMAP_FORMAT) --mask_annotation_files $(CROSSMAP)

gen_count_vectors : gen_p_offsets.txt transcripts_chosen.bed
	get_count_vectors $@ --count_files gen_reads.bam \
		--fiveprime_variable --offset gen_p_offsets.txt \
		--annotation_format BED --annotation_files transcripts_chosen.bed \
		--mask_annotation_format $(CROSSMAP_FORMAT) --mask_annotation_files $(CROSSMAP)

gen_reads_fiveprime_variable_fw.wig gen_reads_fiveprime_variable_rc.wig : gen_reads.bam gen_p_offsets.txt
	make_wiggle --count_files gen_reads.bam \
				--fiveprime --offset 0 \
				-o gen_reads_fiveprime_0
	make_wiggle --count_files gen_reads.bam \
				--threeprime --offset 15 \
				-o gen_reads_threeprime_15
	make_wiggle --count_files gen_reads.bam \
				--threeprime --offset 0 \
				-o gen_reads_threeprime_0
	make_wiggle --count_files gen_reads.bam \
				--center --nibble 0 \
				-o gen_reads_center_0
	make_wiggle --count_files gen_reads.bam \
				--center --nibble 12 \
				-o gen_reads_center_12
	make_wiggle --count_files gen_reads.bam \
				--fiveprime_variable --offset gen_p_offsets.txt \
				-o gen_reads_fiveprime_variable
	make_wiggle --count_files gen_reads.bam \
				--fiveprime --offset 15 \
				-o gen_reads_fiveprime_15

gen_phasing.txt : gen_reads.bam gen_p_offsets.txt transcripts_chosen.bed
	phase_by_size gen --annotation_files transcripts_chosen.bed \
				  --annotation_format BED \
				  --min_length 26 --max_length 31 \
				  --count_files gen_reads.bam \
				  --fiveprime_variable --offset gen_p_offsets.txt
