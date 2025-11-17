#!/bin/bash
#SBATCH --job-name=WH_PED_ASHKENAZI_GOI
#SBATCH --output=/hpc/umc_laat/gvandersluis/outs_sbatch/%x_%j.out
#SBATCH --error=/hpc/umc_laat/gvandersluis/errors_sbatch/%x_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
# Email notifications
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=g.vandersluis-2@umcutrecht.nl


apptainer exec \
 --bind /hpc/umc_laat \
 /hpc/umc_laat/gvandersluis/software/whatshap_v1.sif \
 whatshap phase \
 --reference=/hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
 --ped=/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/ashk_giab.ped \
 -o /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/phased_ashk_ROI.bcf \
 /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/trio_ashk_ROI.bcf \
 /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/calls.sorted.bam \
 /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG003/calls.sorted.bam \
 /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG004/calls.sorted.bam
