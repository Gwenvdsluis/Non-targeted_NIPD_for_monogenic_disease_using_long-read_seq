#!/bin/bash
#SBATCH --job-name=WF_HUMAN_VARIATION_HG2
#SBATCH --output=/hpc/umc_laat/gvandersluis/outs_sbatch/%x_%j.out
#SBATCH --error=/hpc/umc_laat/gvandersluis/errors_sbatch/%x_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
# Email notifications
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=g.vandersluis-2@umcutrecht.nl

# --- PATH fix so `singularity` exists ---
export PATH="$HOME/bin:$PATH"

# --- Initialise SDKMAN ---
export SDKMAN_DIR="$HOME/.sdkman"
source "$SDKMAN_DIR/bin/sdkman-init.sh"

# --- Force Java 17 (Temurin) ---
sdk use java 17.0.17-tem

export NXF_SINGULARITY_ENGINE=apptainer
export NXF_SINGULARITY_CACHEDIR="$HOME/.nextflow/singularity"

#export NXF_SINGULARITY_ENGINE=apptainer

data_file="/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/calls.sorted.bam"
output_dir="/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/WF_H_VAR/output"
work_dir="/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/WF_H_VAR"
ref_dir="/hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
hwf_sing="/hpc/umc_laat/gvandersluis/software/wf-human-variation_latest.sif"


NXF_VER=23.10.1 nextflow run epi2me-labs/wf-human-variation \
	--bam $data_file \
	--ref $ref_dir \
	--sample_name 'paw70337_try' \
	--snp \
	--sv \
	--mod \
	--phased \
        --igv \
	-profile singularity \
        --out_dir $output_dir 
