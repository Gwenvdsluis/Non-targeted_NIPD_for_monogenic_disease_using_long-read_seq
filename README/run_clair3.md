### README ###

Testing the variant calling tool CLAIR3

BASIC INPUTS to CLAIR3
- BAM file (basecalling ONT)
- Refgenome

--bam_fn :
From "https://42basepairs.com/browse/s3/ont-open-data/giab_2025.01/basecalling/sup"
the HG002 (son), HG003 (father) and HG004 (mother) data was used from the Ashkenazi trio (Only the top sample).

For Clair3, a .bam file for each of the individuals was downloaded to a folder dedicated to their sample:
- data/Ont_data_nhung/HG002
- data/Ont_data_nhung/HG003
- data/Ont_data_nhung/HG004

--ref_fn :
The reference genome was donwloaded from "https://42basepairs.com/browse/s3/ont-open-data/gm24385_2020.09/config/ref"
this is saved in /hpc/umc_laat/resources/refgenomes/

--threads :
Depends on own machine settings, I used 4

--platform :
I used "ont" as ONT data is used

--model_path :
CLAIR3 model path

--output : 
the output is saved to folder inside of the sample folders:
- data/Ont_data_nhung/HG002/variant_calls_HG002/
- data/Ont_data_nhung/HG003/variant_calls_HG003/
- data/Ont_data_nhung/HG004/variant_calls_HG004/

EXTRA INPPUT CLAIR3

--bed_fn : 
chromosomal locations of interested area in a bed file
This file was made via filter_hg.sh
This script uses the annotated human genome from gencode and filters out only the interested genes (HBB and CTFR)
