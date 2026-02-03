#!/bin/bash
#SBATCH --job-name=HAPLOTAG_4
#SBATCH --output=/hpc/umc_laat/gvandersluis/outs_sbatch/%x_%j.out
#SBATCH --error=/hpc/umc_laat/gvandersluis/errors_sbatch/%x_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=8
# Email notifications
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=g.vandersluis-2@umcutrecht.nl

samp="4"
vers="wf_h_var"

apptainer exec --bind /hpc/umc_laat/ \
  /hpc/umc_laat/gvandersluis/software/whatshap_v1.sif \
  whatshap haplotag -o /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00${samp}/${vers}_ROI/haplotagged.bam --reference /hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --ignore-read-groups \
  /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00${samp}/${vers}/${vers}.wf_snp.vcf.gz /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00${samp}/calls.sorted.bam
