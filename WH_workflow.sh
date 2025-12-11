#!/bin/bashsx
#SBATCH --job-name=WH_workflow_HG2_OMIM
#SBATCH --output=/hpc/umc_laat/gvandersluis/outs_sbatch/%x_%j.out
#SBATCH --error=/hpc/umc_laat/gvandersluis/errors_sbatch/%x_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
# Email notifictions
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=g.vandersluis-2@umcutrecht.nl


sample="2"
version="OMIM"
folder="/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00${sample}"
vcf="${folder}/SAMPLE.wf_snp.vcf.gz"
bam="${folder}/calls.sorted.bam"
ref_G="/hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
roi="/hpc/umc_laat/gvandersluis/data/Ref_HG/HG_OMIM_ROI.bed"
phased_vcf="${folder}/${version}_ROI/phased_ROI_HG${sample}.vcf"
benchmark="${folder}/HG00${sample}_BM_SSANDT_rn.vcf.gz"

TMPDIR="/hpc/umc_laat/gvandersluis/tmp"
mkdir -p "${TMPDIR}"

mkdir -p "${folder}/${version}_ROI"

# To match the sample names of the vcf to the bam files, the following command was used for each sample:
echo "=== Reheader VCF ==="
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bcftools reheader -s <( echo hg00${sample} ) -o ${folder}/SAMPLE_renamed.vcf.gz ${vcf}

echo " index vcf "
# Index
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bcftools index -t ${folder}/SAMPLE_renamed.vcf.gz


# Filter the vcf file for QUAL > 40 and the defined regions of interest of the OMIM genes
echo "=== FIlter VCF for QUAL >= & ROI ==="
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bash -c "
    bcftools view -f PASS -i 'QUAL>=40' -R ${roi} ${folder}/SAMPLE_renamed.vcf.gz | bcftools sort -T ${TMPDIR} -Oz -o ${folder}/${version}_ROI/Sample_${version}_ROI_d.vcf.gz"

echo " DROP DUPLICATES "
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bcftools norm -D ${folder}/${version}_ROI/Sample_${version}_ROI_d.vcf.gz -o ${folder}/${version}_ROI/Sample_${version}_ROI.vcf.gz

# Index
echo " index vcf "
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bcftools index -t ${folder}/${version}_ROI/Sample_${version}_ROI.vcf.gz

# RUN WHATSHAP PHASE
  # inputs: reference genome, filtered vcf, bam file
echo "=== RUN WHATSHAP PHASE ==="
apptainer exec \
    --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/whatshap_v1.sif \
    whatshap phase -o ${phased_vcf} \
    --reference=${ref_G} ${folder}/${version}_ROI/Sample_${version}_ROI.vcf.gz ${bam}

echo " make a no header version of the vcf for analysis "
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bcftools view -H ${phased_vcf} > ${folder}/${version}_ROI/phased_ROI_HG${sample}_nh.vcf


# Index the phased vcf file
echo "ZIP phased vcf" 
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bcftools view -Oz -o ${phased_vcf}.gz ${phased_vcf}

echo " index phased vcf"
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bcftools index -t ${phased_vcf}.gz

# Make GTF file
echo "make GTF file"
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/whatshap_v1.sif \
    whatshap stats --gtf=${folder}/${version}_ROI/phased_${version}_ROI_hg${sample}.gtf ${phased_vcf}

# Make stats TSV
echo "Make TSV file"
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/whatshap_v1.sif \
    whatshap stats --tsv=${folder}/${version}_ROI/stats_phased_${version}_ROI_hg${sample}.tsv ${phased_vcf}

echo "=== Haplotag the BAM ==="
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/whatshap_v1.sif \
    whatshap haplotag -o ${folder}/${version}_ROI/haplotagged.bam --reference ${ref_G} ${phased_vcf}.gz ${bam}

# Prep the BM file
echo "=== Reheader the Benchmark file ==="
apptainer exec --bind /hpc/umc_laat/ \
    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
    bcftools reheader -s <(echo hg00${sample}) -o ${folder}/HG00${sample}_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz ${benchmark}

#echo " index benchmark"
#apptainer exec --bind /hpc/umc_laat/ \
#    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
#    bcftools index -t ${benchmark}

# SWITCH ERRORS
#apptainer exec --bind /hpc/umc_laat/ \
#    /hpc/umc_laat/gvandersluis/software/whatshap_v1.sif \
#    whatshap compare --tsv=eval_${version)_HG${sample}.tsv $ ${phased_vcf}

# no header version of the phased VCF
#apptainer exec --bind /hpc/umc_laat/ \
#    /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif \
#    bcftools view -H ${phased_vcf} ${folder}/${version}_ROI/phased_${version}_ROI_hg${sample}_nh.vcf
#
