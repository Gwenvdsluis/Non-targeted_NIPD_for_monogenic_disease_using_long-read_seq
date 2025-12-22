data_file="/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/calls.sorted.bam"
output_dir="/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/WF_H_VAR/output"
work_dir="/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/WF_H_VAR"
ref_dir="/hpc/umc_laat/resources/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"



nextflow run epi2me-labs/wf-human-variation \
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
