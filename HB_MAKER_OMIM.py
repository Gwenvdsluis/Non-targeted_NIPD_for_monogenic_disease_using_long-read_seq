#!/usr/bin/env python3

### Imports

import pandas as pd
import numpy as np
import subprocess
import os
import pysam


def create_directories(sample, version, tmpdir):
    ### Ensure TMPDIR exists
    os.makedirs(tmpdir, exist_ok=True)
    ### Ensure the output file exists
    out_dir = f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{sample}/{version}_ROI"
    os.makedirs(out_dir, exist_ok=True)


def ROI_getter(phased_vars, roi_b, tmpdr, phased_vars_roi, version, version_bm, sample):
    print("ROI_GETTER")
    ### Get the ROI of the OMIM genes in the phased vcf
    cmd_roi = f"apptainer exec -B /hpc/:/hpc/ /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif bash -c 'bcftools view -R {roi_b} {phased_vars} | bcftools sort -T {tmpdr} -Oz -o {phased_vars_roi}'"
    ### Run the command
    subprocess.run(cmd_roi, shell=True, check=True)
    ### Index the ROI VCF
    cmd_index = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif", "bcftools", "index", "-t", phased_vars_roi]
    subprocess.run(cmd_index, check=True)
    print("DONE")
    print("Reheader:")
    cmd_bm_reheader = f"apptainer exec -B /hpc/:/hpc/ /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif bcftools reheader -s <(echo {version_bm}) -o /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{sample}/{version}_ROI/HG00{sample}_BM_SSANDT_rn.vcf /hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{sample}/HG00{sample}_GRCh38_1_22_v4.2.1_benchmark_phased_MHCassembly_StrandSeqANDTrio.vcf.gz"
    subprocess.run(cmd_bm_reheader, shell=True, executable="/bin/bash", check=True)
    print("Rheader done")
    bm_index = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif", "bcftools", "index", "-t", f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{sample}/{version}_ROI/HG00{sample}_BM_SSANDT_rn.vcf"]
    subprocess.run(bm_index, check=True)
    print("Index done")


def hb_maker(phased_vars_file):
    ### Read the filtered vcf file (filtered on ROI of the OMIM genes & QUAL >= 40) and get the desired columns
    df_vcf = pd.read_csv(phased_vars_file, sep="\t", comment='#', names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', "sample"])
    df_vcf["GENE"] = df_vcf['INFO'].str.split("|").str[4]
    df_splitted = df_vcf['sample'].str.split(":")
    df_vcf["GT"] = df_splitted.str[0]
    df_vcf["GQ"] = df_splitted.str[1]
    df_vcf["DP"] = df_splitted.str[2]
    df_vcf["PS"] = df_splitted.str[5]

    ### Only keep the phased variants (the variants with a phase tag)
    phased = df_vcf.loc[df_vcf['PS'].notna(), ['CHROM','POS','QUAL','GENE','GT','GQ','DP','PS','REF','ALT']]
    ### Group by phase tag and get the min, max, and amount of variants per haploblock
    phased['POS'] = pd.to_numeric(phased['POS'], errors='coerce')
    phased['GQ'] = pd.to_numeric(phased['GQ'], errors='coerce')
    grouped = phased.groupby(['PS']).agg({'POS': ['min','max', 'count'],
                                                    'CHROM': 'first',
                                                    'GQ': 'mean',
                                                    'GENE': 'unique'}).reset_index()

    ### Change the column names, type, and add desired columns
    grouped.columns = ['_'.join(col).rstrip('_') if isinstance(col, tuple) else col for col in grouped.columns]
    grouped.rename(columns={"POS_min":"START_HB", "POS_max": "END_HB", "POS_count":"phased_variants", "CHROM_first":"chromosome", "GENE_unique":"GENES_in_HB"}, inplace=True)
    #grouped["PS_tag"] = grouped["PS_tag"].astype(int)
    grouped["HB_length"]=grouped["END_HB"]-grouped["START_HB"]
    grouped["total_VARS"] = None
    grouped["PS_tag"] = grouped["chromosome"]+"_"+grouped["PS"]
    print(grouped)
    return grouped


def run_switch(haploblock_table, sample, version, ph_vars, total_vars):
    ### A counter to check how many haploblocks are already compared to BM
    cnt = 0 

    ### Variables and Files for the output
    folder = "/hpc/umc_laat/gvandersluis/data/"
    BM_vars = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/HG00{sample}_BM_SSANDT_rn.vcf"
    BM_ROI = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/ROI_eval.tsv"
    open(BM_ROI, "w").close()
    switches_ROI = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/switches_ROI.bed"
    open(switches_ROI, "w").close()

    ### Loop through every haploblock
    for idx, item in haploblock_table.iterrows():
        ### Get the haploblock region
        reg = item["chromosome"]+":"+str(item["START_HB"])+"-"+str(item["END_HB"])
        print(item["PS_tag"])
        print(reg)
        ### ADD TOTAL VARIANTS in the region from the unfiltered vcf file
        cmd = f"apptainer exec -B /hpc/:/hpc/ /hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif bcftools view -H -r {reg} {total_vars} | wc -l"
        ### Run the command
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        ### Get the count as integer and add the value to the column
        count = int(result.stdout.strip())
        haploblock_table.loc[idx, "total_VARS"] = count

        ### WHATSHAP COMPARE
        PStag = str(item["PS_tag"])
        if "NaN" not in reg:
            ### Prepare the benchmark file for whatshap compare
            ### make an output file for the filtered Benchmark file (filtered on region of interest) and run that code + index the vcf
            out_vcf = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/{PStag}_BM.vcf"
            cmd = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif", "bcftools", "view", "-r", reg, BM_vars, "-Oz", "-o", out_vcf]
            subprocess.run(cmd, check=True)
            indexer = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif", "bcftools", "index", "-t", out_vcf]
            subprocess.run(indexer, check=True)

            ### Make specific region version of the phased vcf file (Prepare file for whatshap compare)
            out_ph_vars_vcf = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/{PStag}_phased_vars.vcf"
            cmd_1 = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif", "bcftools", "view", "-r", reg, ph_vars, "-Oz", "-o", out_ph_vars_vcf]
            r = subprocess.run(cmd_1, check=True)
            print(r.stderr, "\n", r.stdout)
            indexer_1 = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif", "bcftools", "index", "-t", out_ph_vars_vcf]
            subprocess.run(indexer_1, check=True)

            ### Prepare output files for WHATSHAP compare
            indexed = out_vcf+".tbi"
            indexed_1 = out_ph_vars_vcf+".tbi"
            tsv_out = f"{PStag}_eval.tsv"
            open(tsv_out, "a+").write("")

            ### If the output vcf file has variants (If there are variants found in the benchmark file of this particular region, then:)
            if sum(1 for _ in pysam.VariantFile(out_vcf)) != 0:
                ### Make file for the switch error bed file
                bed_out = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/{PStag}_switch.bed"
                ### RUN WHATSHAP COMPARE
                cmd2 = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/whatshap_v1.sif", "whatshap", "compare","--switch-error-bed", bed_out, "--tsv-pairwise", tsv_out, "--names", f"BENCHMARK,{PStag}", out_vcf, out_ph_vars_vcf]
                subprocess.run(cmd2, check=True, capture_output=True)

                ### If the whatsap compare output file is empty (first run), then write output
                if os.path.getsize(BM_ROI) == 0:
                    with open(BM_ROI, "w") as out:
                        with open(tsv_out, "r") as inp:
                            out.write(inp.read())
                ### If whatshap compare output file is not empty anymore (every run after the first), then append output
                else:
                    with open(BM_ROI, "a") as out:
                        with open(tsv_out, "r") as inp:
                            next(inp) # Skip header
                            out.write(inp.read())
                ### Write or append switch error locations to the bed file
                with open(switches_ROI, "a+") as out_b:
                    with open(bed_out, "r") as inp_b:
                        out_b.write(inp_b.read())
            ### REMOVE temporary needed files
                os.remove(bed_out)
            os.remove(out_vcf)
            os.remove(indexed)
            os.remove(out_ph_vars_vcf)
            os.remove(indexed_1)
            os.remove(tsv_out)
        print(PStag, "\t Done")

        cnt += 1
        print(cnt)
    return BM_ROI, haploblock_table

def format_haploblock_table(bench_roi, hb_comp, sam, ver):
    sw_e = pd.read_csv(bench_roi, sep="\t")

    ### Format the switch error file
    sw_e = sw_e[["dataset_name1","het_variants0","all_switches","all_switch_rate","all_switchflips","all_switchflip_rate","blockwise_hamming_rate"]]
    sw_e["Accuracy"] = (1-sw_e["blockwise_hamming_rate"]).astype(str)+"%"
    sw_e.drop(columns="blockwise_hamming_rate")
    sw_e.rename(columns={"dataset_name1": "PS_tag"}, inplace=True)
    sw_e = pd.merge(hb_comp, sw_e, on="PS_tag", how="left")
    sw_e.sort_values("all_switches").tail(6)

    sw_e.to_csv(f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{sam}/{ver}_ROI/Haploblock_switches.csv", index=False)

### RUN scripts
def main():
    ### RUNNING ON NEW DATA:
    # change: 'vers', 'samp', & filenames (phased_variants & phased_variants_roi)

    print("HI, starting")
    samp="2"
    vers="SUP_v5.2"
    reheader_bm_s_name="SUP_v5.2"
    TMPDIR="/hpc/umc_laat/gvandersluis/tmp"

    ### Create necessary directories
    create_directories(samp, vers, TMPDIR)

    ### Files needed for the analysis
    phased_variants=f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{samp}/{vers}/SUP_v5.2.wf_snp.vcf.gz" # if filter; dont change this file
    roi_bed=f"/hpc/umc_laat/gvandersluis/data/Ref_HG/HG_OMIM_ROI_merged.bed"
    phased_variants_roi=f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{samp}/{vers}_ROI/SUP_v5.2.wf_snp.vcf.gz" # if filter; change this file

    ### Create phased vcf of only the ROI
    ROI_getter(phased_variants, roi_bed, TMPDIR, phased_variants_roi, vers, reheader_bm_s_name, samp)

    ### Generate the VCF haploblock table of the ROI
    grouped = hb_maker(phased_variants_roi)
    print("First function DONE")

    ### Whatshap compare per haploblock
    bm_roi, hb_tab = run_switch(grouped, samp, vers, phased_variants_roi, phased_variants)
    hb_tab.to_csv(f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG00{samp}/{vers}_ROI/final_df.csv", index=False)

    format_haploblock_table(bm_roi, hb_tab, samp, vers)


main()
