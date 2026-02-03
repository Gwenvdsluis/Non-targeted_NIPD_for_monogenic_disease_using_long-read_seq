#!/bin/bash
#SBATCH --job-name=check_cor_phased_HG3
#SBATCH --output=/hpc/umc_laat/gvandersluis/outs_sbatch/%x_%j.out
#SBATCH --error=/hpc/umc_laat/gvandersluis/errors_sbatch/%x_%j.err
#SBATCH --time=5:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
# Email notifications
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=g.vandersluis-2@umcutrecht.nl


python3 -u /hpc/umc_laat/gvandersluis/scripts/Corr_phased.py
