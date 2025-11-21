#!/usr/bin/env python3

import subprocess
import os

sample = "2"
version = "filtered"
folder = "/hpc/umc_laat/gvandersluis/data/"
ph_vars = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/phased_ROI.vcf.gz"
roi = f"{folder}Ref_HG/HG_annotation_ROI.bed"

BM_ROI = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/ROI_eval.tsv"
open(BM_ROI, "w").close()


for region in open(roi, "r").readlines():
    regspl = region.split("\t")
    gen = regspl[3]
    reg = regspl[0]+":"+regspl[1]+"-"+regspl[2]
    out_vcf = f"{folder}Ont_data_nhung/HG00{sample}/{version}_ROI/{gen}_BM.vcf"
    cmd = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/bcftools_v1.9-1-deb_cv1.sif", "bcftools", "view", "-r", reg, ph_vars, "-Oz", "-o", out_vcf]
    subprocess.run(cmd, check=True)
    tsv_out = f"{gen}_eval.tsv"
    cmd2 = ["apptainer", "exec", "-B", "/hpc/:/hpc/", "/hpc/umc_laat/gvandersluis/software/whatshap_v1.sif", "whatshap", "compare", "--tsv-pairwise", tsv_out, "--names", f"BENCHMARK,{gen}", out_vcf, ph_vars]
    subprocess.run(cmd2, check=True)
    if os.path.getsize(BM_ROI) == 0:
        with open(BM_ROI, "w") as out:
            with open(tsv_out, "r") as inp:
                out.write(inp.read())
    else:
        with open(BM_ROI, "a") as out:
            with open(tsv_out, "r") as inp:
                next(inp)
                out.write(inp.read())
    os.remove(out_vcf)
    os.remove(tsv_out)
    print(gen, "\t Done")

