#!/usr/bin/env python3

# loop through the haploblock file
# for every row make line for the newly generated GTF file
HBs_f= "/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/wf_hvar_ROI/Haploblock_switches.csv"
HBs = open(HBs_f, "r").readlines()
out_gtf =open("/hpc/umc_laat/gvandersluis/data/Ont_data_nhung/HG002/wf_hvar_ROI/phased_vars_hblks.gtf", "w")

for hb in HBs[1::]:
    new_line=hb.split(",")[4]+'\tPhasing\texon\t'+hb.split(",")[1]+'\t'+hb.split(",")[2]+'\t.\t+\t.\tgene_id "'+hb.split(",")[0]+'"; transcript_id "'+hb.split(",")[0]+'.1";\n'
    out_gtf.write(new_line)

out_gtf.close()
