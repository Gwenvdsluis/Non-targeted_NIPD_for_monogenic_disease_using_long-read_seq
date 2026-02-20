### README ###

### Performing the whole pipeline from scratch

## Basecalling file from Epi2me

## Run the ONT WF :
# you need the:
# calls.sorted.bam # BAM FILE
# WF_H_VAR # Directory in which it can be saved
# Reference genome
# Sample name and email-adress
 
# sh run.sh /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/calls.sorted.bam \
# /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/WF_H_VAR/ \
# /hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
# wf_hvar g.vandersluis-2@umcutrecht.nl

## run_HB_MAKER_OMIM.sh
# This script makes the HB table (including the sub-haploblocks) and performs whatshap compare

## run_corr_phased.sh
# This script makes a file that defines the variants to be called incorrectly, correctly or unknown

## Haplotagging
# via run_realign_sort.sh if the bam file is not aligned yet
# via run_hapltotag.sh if it is

## FOR visualisations
# filter variants for only phased variants that are also phased in the benchmark & whithin the ROI
# intersect with bcftools to identify all variants on one haplotype
# bcftools isec  -p dir roi_phased.vcf.gz phased_BM.vcf.gz
