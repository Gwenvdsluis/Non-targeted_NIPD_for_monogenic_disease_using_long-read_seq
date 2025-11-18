### README ###

Run the WhatsHap tool on pedigree mode, using pedigree information for mapping.

BASIC INPUTS
- reference genome (.fasta)
- pedigree file (.ped)
- merged trio (.bcf)
- base call files (.bam) 


THE REFERENCE GENOME

The reference genome was donwloaded from "https://42basepairs.com/browse/s3/ont-open-data/gm24385_2020.09/config/ref"
this is saved in /hpc/umc_laat/resources/refgenomes/


PEDIGREE FILE

self generated .ped file with a family identifier. This file defines the relations between the family members,
including their gender and if they are affected (optional, not used here)
The file is saved in 
data/Ont_data_nhung/ashk_giab.ped


MERGED TRIO (merging the vcf files of the ped samples)

From "https://42basepairs.com/browse/s3/ont-open-data/giab_2025.01/analysis/wf-human-variation/sup"
the HG002 (son), HG003 (father) and HG004 (mother) vcf files were downloaded from the Ashkenazi tri$

A .vcf.gz file for each of the individuals was downloaded to a folder dedicated to their sample:
- data/Ont_data_nhung/HG002/SAMPLE_.wf_snp.vcf.gz
- data/Ont_data_nhung/HG003/SAMPLE_.wf_snp.vcf.gz
- data/Ont_data_nhung/HG004/SAMPLE_.wf_snp.vcf.gz

To match the sample names of the vcf to the bam files, the following command was used for each sample:
bcftools reheader -s <(echo hg004) -o data/Ont_data_nhung/HG004/SAMPLE_renamed.vcf.gz data/Ont_data_nhung/HG004/SAMPLE.wf_snp.vcf.gz

Resulting in a SAMPLE_renamed.vcf.gz file in each sample folder

Next the vcf files were filtered for only variants in the specified genes that passed (PASS) the filtering steps:
bcftools view -f PASS -i 'QUAL>=40' -R /data/Refgen/HG_annotated_ROI.bed SAMPLE_renamed.vcf.gz -Oz -o Sample_filtered_ROI.vcf.gz

Resulting in a Sample_filtered_ROI.vcf.gz file for each sample in the designated folder mentioned previously

These files were merged into one .bcf file using:

bcftools merge -o trio_ashk_ROI.bcf HG002/Sample_filtered_ROI.vcf.gz HG003/Sample_filtered_ROI.vcf.gz HG004/Sample_filtered_ROI.vcf.gz

bcftools index trio_ashk_ROI.bcf

BASECALLING FILES

From "https://42basepairs.com/browse/s3/ont-open-data/giab_2025.01/basecalling/sup"
the HG002 (son), HG003 (father) and HG004 (mother) bam files were downloaded from the Ashkenazi trio (Only the top sample).

A .bam file for each of the individuals were downloaded to a folder dedicated to their sample:
- data/Ont_data_nhung/HG002
- data/Ont_data_nhung/HG003
- data/Ont_data_nhung/HG004

via run_ped_whatshap.sh:
The output can be found in: data/Ont_data_nhung/phased_ashk_ROI.bcf

INDEX 
	bcftools index phased_ashk_ROI.bcf

MAKE GTF
	whatshap stats --gtf=phased_ashk_ROI.gtf phased_ashk_ROI.bcf

MAKE STATS TSV
	whatshap stats --tsv=stats_phased_ROI_trio.tsv phased_ashk_ROI.bcf
