#!/bin/bash

### Gwen van der Sluis
### Bash 4.4.20(1)-release
### GNU nano 2.9.8
### 27-10-25

####
# This script filteres the Genes Of Interest (GOI) from the hg annotation file from gencode.
# First only coding genes are filterd out,
# Then the GOI are selected
# At last the locations are selected with a start and end position 1Mbp further
# if the start position is a negative value it is set to 1, if it is further than the max chromosome lenght it is set to the
# max chromosome lenght
####

HG="/hpc/umc_laat/gvandersluis/data/Ref_HG"

# Chromosome length dictionary
declare -A chr_l=(["chr1"]="248956422" ["chr2"]="242193529" ["chr3"]="198295559" ["chr4"]="190214555" ["chr5"]="181538259" \
["chr6"]="170805979" ["chr7"]="159345973" ["chr8"]="145138636" ["chr9"]="138394717" ["chr10"]="133797422" \
["chr11"]="135086622" ["chr12"]="133275309" ["chr13"]="114364328" ["chr14"]="107043718" ["chr15"]="101991189" \
["chr16"]="90338345" ["chr17"]="83257441" ["chr18"]="80373285" ["chr19"]="58617616" ["chr20"]="64444167" \
["chr21"]="46709983" ["chr22"]="50818468" ["chrX"]="156040895" ["chrY"]="57227415")

# Function that loops through the dictionary
awk_dict_init=$(for k in "${!chr_l[@]}"; do
  printf 'chr_l["%s"]=%s; ' "$k" "${chr_l[$k]}"
done)



zcat $HG/gencode.v49.annotation.gtf.gz | \
 awk '($3=="gene")' | \
 grep -E '"HBB";|"BRCA2";|"BRCA1";|"POLG";|"TSEN54";|"CRTAP";|"SMN1";|"PEX7";|"CFTR";|"MUSK";|"CYP21A2";|"HBA1";' | \
 awk -F '\t' '{ if (match($9, /gene_name "([^"]+)"/, a)) print $1"\t"$4"\t"$5"\t"a[1]"\t"$6"\t"$7}' > $HG/HG_annotation_GOI.bed



cat $HG/HG_annotation_GOI.bed | \
  awk -F '\t' 'BEGIN{'"$awk_dict_init"'} # Use the dictionary funtion
   { middle=(($3-$2)/2) + $2
   start=middle-500000; end=middle+500000 # If it is on the plus strand, calculate new start and end position 
   if (start < 1) start=1; # If start position is lower than 1, set it to 1
   if (end > chr_l[$1]) end=chr_l[$1]; # Uses the dict funciton: If value is over the max chormosome length, set to max chr length
   printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1, start, end, $4, $5, $6; # Select the columns neccessary for a bed file
   }' > $HG/HG_annotation_ROI.bed # save to designated file

