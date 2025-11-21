#!/usr/bin/env python3

import pandas as pd
import numpy as np

def hb_maker(phased_vars_file):
    df_vcf = pd.read_csv(phased_vars_file, sep="\t", names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', "sample"])
    df_vcf["GENE"] = df_vcf['INFO'].str.split("|").str[4]
    df_splitted = df_vcf['sample'].str.split(":")
    df_vcf["GT"] = df_splitted.str[0]
    df_vcf["DP"] = df_splitted.str[2]
    df_vcf["PS"] = df_splitted.str[5]
 
    # Min Quality of 10, this means there is no more than 1/10000 chance of an icorrect base call (>99.9% accuracy)
    # & (pd.to_numeric(df_vcf['DP']) >= 20
    phased = df_vcf.loc[df_vcf['PS'].notna(), ['CHROM','POS','QUAL','GENE','GT','DP','PS','REF','ALT']]
    phased["DP"] = pd.to_numeric(phased["DP"], errors="coerce")
    phased = phased.sort_values(by="DP")

    #print(phased[phased["GENE"].str.contains("SM")])


    gr_unfilterd = df_vcf.groupby(['PS']).agg({'POS': ['min','max', 'count'],
                                                    'CHROM': 'first',
                                                    'QUAL': 'mean',
                                                    'GENE': 'unique'}).reset_index()

    grouped = phased.groupby(['PS']).agg({'POS': ['min','max', 'count'],
                                                    'CHROM': 'first',
                                                    'QUAL': 'mean',
                                                    'GENE': 'unique'}).reset_index()
    #print(grouped[grouped[("CHROM", "first")]=="chr7"])

    #print(ph_GOI[ph_GOI[("POS")]==117524809])
    #print(gr_unfilterd)
    #print(grouped)
    #mr_grouped = pd.merge(grouped, gr_unfilterd, how="left", on=["PS",("POS","max")])
    print(len(phased))
    return grouped


def phased_GOI(regions_file, phased_hps):
    genes_of_interest = pd.read_csv(regions_file, sep="\t", names=["CHROM","START","END","GENE","INFO","STRAND"])
    goi = genes_of_interest["GENE"].unique()
    genes_of_interest["POSITION"] = genes_of_interest["CHROM"].astype(str)+":"+genes_of_interest["START"].astype(str)+"-"+genes_of_interest["END"].astype(str)
    #print(genes_of_interest)
    # Make sure each GENE entry is a Python list
    phased_hps["GENE"] = phased_hps[("GENE","unique")].apply(lambda x: x.tolist())
    #print(phased_hps)

    # Filter rows where at least one gene in 'unique' is in GOI
    filtered = phased_hps[phased_hps[("GENE", "unique")].apply(lambda genes: any(g in goi for g in genes))].copy()

    # Now, inside each row, keep only the genes that match the GOI
    filtered[("GENE", "unique")] = filtered[("GENE", "unique")].apply(lambda genes: ", ".join([str(g) for g in genes if g in goi]))
    # Remove sub columns
    filtered.columns = ['_'.join([c for c in col if c]).strip() if isinstance(col, tuple) else col for col in filtered.columns.values]
    # Now rename for clarity
    filtered = filtered.rename(columns={'PS': 'PS',
                                        'POS_min': 'POS_START',
                                        'POS_max': 'POS_END',
                                        'POS_count': 'PHASED_VARIANTS',
                                        'CHROM_first': 'chromosome',
                                        'QUAL_mean': 'QUAL_MEAN',
                                        'GENE_unique': 'GENE'})

    filtered["HB_LENGTH"] = filtered["POS_END"]-filtered["POS_START"]
    # print(filtered)
    # print(genes_of_interest)
    merged = pd.merge(filtered, genes_of_interest, on="GENE", how="right")[["GENE","POSITION","POS_START","POS_END","HB_LENGTH","PHASED_VARIANTS"]]

    return merged




    # Now select only the columns you want to display
    # merged = merged[
    #     [
    #     "GENE",
    #     "CHROM_phased",  # from filtered
    #     "POS_min", "POS_max", "QUAL_mean",
    #     "CHROM_goi", "START", "END", "STRAND"  # from genes_of_interest
    #     ]
    # ]
    # print(merged)
    # GOI_variants = pd.merge(phased_hps, genes_of_interest, how="inner",on=["CHROM","GENE"])
    # print(GOI_variants)
    

def switch_errors(switch_err_file, bs_df):
    sw_e = pd.read_csv(switch_err_file, sep="\t")
    sw_e= sw_e[["dataset_name0","all_switches"]]
    sw_e = sw_e.rename(columns={"dataset_name0": "GENE"})

    switch_df = pd.merge(bs_df, sw_e, on="GENE", how="inner")
    print(switch_df, '\n\n')
    return switch_df




def main():
    goi_file = "/hpc/umc_laat/gvandersluis/data/Ref_HG/HG_annotation_ROI.bed"
    filtered_v = ["/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/filtered_ROI/phased_ROI_nh.vcf",
                "/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/filtered_ROI/ROI_eval.tsv"]
    
    print("VCF filtered on PASS, QUAL >= 40 & ROI")
    phased_hb = hb_maker(filtered_v[0])
    basis_df = phased_GOI(goi_file, phased_hb)
    sw_er_db1 = switch_errors(filtered_v[1], basis_df)

    print("VCF filtered on PASS & ROI")
    unfiltered_v = ["/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/unfiltered_ROI/uf_phased_ROI_nh.vcf",
                "/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/unfiltered_ROI/ROI_eval.tsv"]
    phased_hb2 = hb_maker(unfiltered_v[0])
    basis_df2 = phased_GOI(goi_file, phased_hb2)
    sw_er_db2 = switch_errors(unfiltered_v[1], basis_df2)

    print("VCF filtered on ROI")
    ununfiltered_v = ["/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/un_unfiltered_ROI/uf_uf_phased_ROI_nh.vcf",
                "/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/un_unfiltered_ROI/ROI_eval.tsv"]
    phased_hb3 = hb_maker(ununfiltered_v[0])
    basis_df3 = phased_GOI(goi_file, phased_hb3)
    sw_er_db3 = switch_errors(ununfiltered_v[1], basis_df3)
    sw_er_db3 = sw_er_db3[["GENE","PHASED_VARIANTS"]]
    sw_er_db3.rename(columns={"PHASED_VARIANTS": "total_vars"}, inplace=True)


    final_db = pd.merge(sw_er_db1, sw_er_db3, on="GENE", how="inner")
    final_db["FILTERED_VARS"] = final_db["total_vars"] - final_db["PHASED_VARIANTS"]
    final_db.pop("total_vars")

    int_cols = ["POS_START", "POS_END", "HB_LENGTH", "PHASED_VARIANTS", "FILTERED_VARS"]

    for c in int_cols:
        final_db[c] = pd.to_numeric(final_db[c], errors="coerce").astype("Int64")
    final_db = final_db.astype(str).replace("<NA>", "NaN")
    print(final_db.dtypes)

    print(final_db[["GENE", "POSITION", "POS_START", "POS_END", "HB_LENGTH", "PHASED_VARIANTS", "FILTERED_VARS","all_switches"]])



main()
