#!/usr/bin/env python3

import pandas as pd
import numpy as np
from cyvcf2 import VCF, Writer


def mod_bm(bmf):
    ### Get the Benchmark file
    truth = pd.read_csv(bmf, comment="#", sep="\t", compression="gzip", names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample"])
    ### Get and generate the necessary columns
    truth["GT_BM"]=truth["sample"].str.split(":").str[0]
    truth["tag"]=truth["CHROM"].astype(str)+"_"+truth["POS"].astype(str)
    short_truth=truth[["tag", "GT_BM"]]
    ### Remove the variants that contain more variation possibilities '2' and / (unphased variants)
    short_truth=short_truth[(short_truth["GT_BM"].str.contains("2") == False) & (short_truth["GT_BM"].str.contains("/") == False)]
    return short_truth


def mod_vcf(vcf_f, sh_truth):
    ### Load the variant files, collect the necessary columns and remove unphased variants
    vcf = pd.read_csv(vcf_f, comment="#", sep="\t", names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample"])
    vcf["GT"] = vcf["sample"].str.split(":").str[0]
    vcf = vcf[vcf["GT"].str.contains("/") == False]
    vcf["Cor_phased"]="unknown"
    vcf["PS_tag"] = vcf["CHROM"].astype(str)+"_"+vcf["sample"].str.split(":").str[5]
    vcf["tag"] = vcf["CHROM"].astype(str)+"_"+vcf["POS"].astype(str)

    ### Merge the Benchmark file with the VCF file on the generated location tag 'tag'
    ### Turn all variants of the vcf that are in the benchmark to True.
    ### For haploblocks without switches, this will mean that all variants are correctly phased that are found in the benchmark file
    complete = pd.merge(vcf, sh_truth, on="tag", how="left")
    complete.loc[complete["GT_BM"].notna(), "Cor_phased"] = "True"
    return complete


def informative_tags(hb_f):
    ### From the haploblock file, identify the PS_tags that contain switches.
    ### For these haploblocks the variants need to be determined if they are correctly phased
    v5_2_new = pd.read_csv(hb_f)
    tags = v5_2_new[(v5_2_new["all_switches"] > 0) & (v5_2_new["all_switches"].notna())].sort_values("phased_variants")["PS_tag"].unique()
    return tags


def corr_phased_vars(complete, ps_tags):
    count = 0
    ### Run through all phases tags (haploblocks) in which switches are present
    for i in ps_tags:
        count += 1
        print(count)
        ### Filter rows of the dataframe that contain the PS tag and are not NaN
        subset = complete[(complete["PS_tag"] == i) & (complete["GT_BM"].notna())]
        ### For variant in this cluster if it is not empty:
        if not subset.empty:
            for index, item in subset.iterrows():
                ### To check if the variants are phased correctly according to the benchmark an on which allele they are,
                ### an X or O is connected to the variant. And - as a check to see if something else is up
                if item["GT"] == item["GT_BM"]:
                    complete.loc[index, "Cor_phased"] = "X"
                elif item["GT"] == item["GT_BM"][::-1]:
                    complete.loc[index, "Cor_phased"] = "O"
                else:
                    complete.loc[index, "Cor_phased"] = "-"
        ### Per haploblock we can now calculate the Hemming rate.
        tot = complete[complete["PS_tag"] == i].groupby("Cor_phased").agg("count")
        tot_x = tot.loc[["X"],"GT"]
        tot_o = tot.loc[["O"],"GT"]
        hemm_rate_x = tot_x.loc["X"]/(tot_x.loc["X"] + tot_o.loc["O"])
        hemm_rate_o = tot_o.loc["O"]/(tot_x.loc["X"] + tot_o.loc["O"])

        ### If the hemming rate is higher with the X, then O variants are the switches and the other way around.
        ### If the rate is 0.50, it is not possible to know which are phased correctly
        if round(hemm_rate_x,2) == 0.50:
            complete.loc[(complete["PS_tag"] == i) & (complete["GT_BM"].notna()), "Cor_phased"] = "UNKNOWN"
        elif round(hemm_rate_x,4) > round(hemm_rate_o,4):
            complete.loc[(complete["PS_tag"] == i) & (complete["Cor_phased"] == "O" ), "Cor_phased"] = "False"
            complete.loc[(complete["PS_tag"] == i) & (complete["Cor_phased"] == "X" ), "Cor_phased"] = "True"
        else:
            complete.loc[(complete["PS_tag"] == i) & (complete["Cor_phased"] == "X" ), "Cor_phased"] = "False"
            complete.loc[(complete["PS_tag"] == i) & (complete["Cor_phased"] == "O" ), "Cor_phased"] = "True"
    complete.to_csv("/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/SUP_v5.2_ROI/cor_phased.csv", sep='\t', index=False)
    return complete


def save_vcf(vcff, outf, comp_tab):
    comp_tab = comp_tab.copy()
    ### Index for fast lookup
    comp_tab["key"] = list(zip(comp_tab.CHROM, comp_tab.POS, comp_tab.REF, comp_tab.ALT))
    lookup = dict(zip(comp_tab.key, comp_tab.Cor_phased))

    ### Open original VCF
    vcf_f = VCF(vcff)
    ### Edit the header to include Cor_phased
    vcf_f.add_format_to_header({
        "ID": "Cor_phased",
        "Description": "indicates if the variant is switch error causing",
        "Type": "String",
        "Number": "1"
    })

    ### Write vcf to output file
    writer = Writer(outf, vcf_f)

    ### For every variant, check dataframe and add Cor_phased Column
    for var in vcf_f:
        key = (var.CHROM, var.POS, var.REF, var.ALT[0])
        value = lookup.get(key, "NA")
        var.set_format("Cor_phased", np.array([value.encode("ascii")], dtype="S"))
        writer.write_record(var)
    writer.close()
    vcf_f.close()

def main():
    vers="HG002"
    ### OMIM:
    # sample="OMIM"
    # vcf_file=f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/{vers}/{sample}_ROI/{sample}.vcf.gz"

    ### WF_H_VAR
    # sample="WF_H_VAR"
    # sm="wf_hvar"

    ### SUPv5.2
    sample="SUP_v5.2"
    sm="SUP_v5.2"

    ### If u calculate OMIM, comment the vcf_file below out!
    bm_file = f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/{vers}/{sample}_ROI/{vers}_BM_SSANDT_rn.vcf"
    vcf_file = f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/{vers}/{sample}_ROI/{sm}.wf_snp.vcf.gz"
    haploblock_file = f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/{vers}/{sample}_ROI/Haploblock_switches.csv"
    output_file = f"/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/{vers}/{sample}_ROI/cor_phased.vcf.gz"
    bench_t = mod_bm(bm_file)
    comp_vcf = mod_vcf(vcf_file, bench_t)
    inf_tags = informative_tags(haploblock_file)
    corr_phased = corr_phased_vars(comp_vcf, inf_tags)
    save_vcf(vcf_file, output_file, corr_phased)


main()
