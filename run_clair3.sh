
#!/bin/bash

#SBATCH --job-name=hg2_variant_calling#1_clair3
#SBATCH --output=/hpc/umc_laat/gvandersluis/outs_sbatch/%x_%j.out
#SBATCH --error=/hpc/umc_laat/gvandersluis/errors_sbatch/%x_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
# Email notifications
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=g.vandersluis-2@umcutrecht.nl

SAMPLE=HG002


# Create TMPDIR if missing:
mkdir -p /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/${SAMPLE}/tmp

apptainer exec \
  --bind /hpc/umc_laat \
  --env TMPDIR=/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/${SAMPLE}/tmp \
  /hpc/umc_laat/gvandersluis/software/clair3_v1.2.0.sif \
  /opt/bin/run_clair3.sh \
  --bam_fn=/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/${SAMPLE}/calls.sorted.bam \
  --ref_fn=/hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --threads=4 \
  --platform="ont" \
  --model_path=/opt/models/r1041_e82_400bps_sup_v500 \
  --output=/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/${SAMPLE}/variant_calls_ROI_${SAMPLE} \
  --bed_fn=/hpc/umc_laat/gvandersluis/data/Ref_HG/HG_annotation_ROI.bed
