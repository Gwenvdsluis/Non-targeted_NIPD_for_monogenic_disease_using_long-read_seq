#!/bin/bash

HG="/hpc/umc_laat/gvandersluis/data/Ref_HG"

awk '($5 != "" && $2 == "gene") {print $5}' \
  $HG/mim2gene.txt \
  | grep -F -f - <( zcat $HG/BIOMART_nc_genes.txt.gz ) | awk '{ if ($6 == 1) {$6="+"} else {$6="-"}; print "chr"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$2}' > $HG/HG_OMIM.bed


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


cat $HG/HG_OMIM.bed | \
  awk -F '\t' 'BEGIN{'"$awk_dict_init"'} # Use the dictionary funtion
   { start=$2; end=$3; middle=((end-start)/2) + start
   if ((end-start) < 100000) {
      start=middle-500000; end=middle+500000 # If it is on the plus strand, calculate new start and end position 
   }
   if (start < 1) start=1; # If start position is lower than 1, set it to 1
   if (end > chr_l[$1]) end=chr_l[$1]; # Uses the dict funciton: If value is over the max chormosome length, set to max$
   printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1, start, end, $4, $5, $6; # Select the columns neccessary for a bed file
   }' > $HG/HG_OMIM_ROI.bed # save to designated file



