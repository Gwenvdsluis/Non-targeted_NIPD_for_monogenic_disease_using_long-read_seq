#!/bin/bash

#SBATCH --job-name=WH_HG2_haplo
#SBATCH --output=/hpc/umc_laat/gvandersluis/outs_sbatch/%x_%j.out
#SBATCH --error=/hpc/umc_laat/gvandersluis/errors_sbatch/%x_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
# Email notifications
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=g.vandersluis-2@umcutrecht.nl

SAMPLE=HG002

apptainer exec \
  --bind /hpc/umc_laat/ \
  /hpc/umc_laat/gvandersluis/software/whatshap_v1.sif \
  whatshap haplotag \
  --reference=/hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --output-haplotag-list=/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/$SAMPLE/table_ROI.tsv \
  -o /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/$SAMPLE/tagged_ROI.bam \
  /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/$SAMPLE/phased_ROI.vcf.gz \
  /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/$SAMPLE/calls.sorted.bam 
