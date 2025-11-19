### README ###

Testing out the WhatsHap phasing tool

BASIC INPUTS

- vcf file with variants
- bam file with basecalling
- ref genome
- bed file with regions of intrest

BED FILE: REGIONS OF INTEREST

from gencode the annotated human genome was downloaded and the regions of interest were put in a bed file via filter_hg.sh
The regions were selected with 1Mbp with the gene in the centre.


THE VARIANT CALLING FILE

From "https://42basepairs.com/browse/s3/ont-open-data/giab_2025.01/analysis/wf-human-variation/sup"
the HG002 (son), HG003 (father) and HG004 (mother) vcf files was downloaded from the Ashkenazi trio

A .vcf.gz file for each of the individuals were downloaded to a folder dedicated to their sample:
- data/Ont_data_nhung/HG002/SAMPLE_.wf_snp.vcf.gz
- data/Ont_data_nhung/HG003/SAMPLE_.wf_snp.vcf.gz
- data/Ont_data_nhung/HG004/SAMPLE_.wf_snp.vcf.gz

To match the sample names of the vcf to the bam files, the following command was used for each sample:
bcftools reheader -s <(echo hg004) -o data/Ont_data_nhung/HG004/SAMPLE_renamed.vcf.gz data/Ont_data_nhung/HG004/SAMPLE.wf_snp.vcf.gz

Resulting in a SAMPLE_renamed.vcf.gz file in each sample folder

Next the vcf files were filtered for only variants in the selected regions -R  that passed (PASS) and had a quality higher than or equal to 40,  the filtering steps:
bcftools view -f PASS -i 'QUAL>=40' -R /data/Ref_HG/HG_annotation_ROI.bed SAMPLE_renamed.vcf.gz -Oz -o Sample_filtered_ROI.vcf.gz

Resulting in a Sample_filtered.vcf.gz file for each sample in the designated folder mentioned previously


THE BASECALLING FILE

From "https://42basepairs.com/browse/s3/ont-open-data/giab_2025.01/basecalling/sup"
the HG002 (son), HG003 (father) and HG004 (mother) bam files were downloaded from the Ashkenazi trio (Only the top sample).

A .bam file for each of the individuals was downloaded to a folder dedicated to their sample:
- data/Ont_data_nhung/HG002
- data/Ont_data_nhung/HG003
- data/Ont_data_nhung/HG004


THE REFERENCE GENOME

The reference genome was donwloaded from "https://42basepairs.com/browse/s3/ont-open-data/gm24385_2020.09/config/ref"
this is saved in /hpc/umc_laat/resources/refgenomes/

Via the run_whatshap.sh file the sampes will be saved

The result is saved in /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/phased_ROI.vcf


INDEX
	index via: bgzip -c phased_ROI.vcf > phased_ROI.vcf.gz
		   tabix -p vcf phased_ROI.vcf.gz

MAKE GTF
	whatshap stats --gtf=phased_ROI_HG2.gtf phased_ROI.vcf

MAKE STATS TSV
	whatshap stats --tsv=stats_phased_ROI_HG2.tsv phased_ROI.vcf	


HAPLOTAGGING
via run_haplotag_wh.sh
	whatshap haplotag --reference=/hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta --output-haplotag-list=table_ROI.tsv -o tagged_ROI.bam phased_ROI.vcf.gz calls.sorted.bam
	
	inputs:
	- refgenome
	- phased vcf 
	- sorted bam file

FIND SWITCH ERRORS

whatshap compare --tsv-pairwise=eval_ROI.tsv ../hg002_BM_sstrio.vcf.gz phased_ROI.vcf



MAKE A NO HEADER VERSION OF THE phased VCF

bcftools view -H phased_ROI.vcf > phased_ROI_nh.vcf
