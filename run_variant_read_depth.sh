#!/bin/bash
#SBATCH --job-name=ALLELE_DEPTH_HG2
#SBATCH --output=/hpc/umc_laat/gvandersluis/outs_sbatch/%x_%j.out
#SBATCH --error=/hpc/umc_laat/gvandersluis/errors_sbatch/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
# Email notifications
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=g.vandersluis-2@umcutrecht.nl


# prevent thread oversubscription
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

python3 -u /hpc/umc_laat/gvandersluis/scripts/variant_read_depth.py
