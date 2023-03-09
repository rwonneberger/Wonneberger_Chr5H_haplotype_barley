#!/bin/bash

#################################################################################################################################
## Extract and filter SNPs in the region of interest on chromosome 5H for all accessions in the IPK genebank barley collection ##
#################################################################################################################################

# This script will require the following programs: Tassel, vcftools, samtools, bcftools, tabix.
# For the publication Tassel 5.2.49 and vcftools 0.1.16 were used

# The 5H region of interest is located at ~70 - 320 Mbp in Barke. 
# This region corresponds to ca. 66.2 - 312.2 Mbp in Morex V3. 
# As a more conservative approach we extracted 67-312 Mbp from Morex V3 and saved to MorexV3_5H_67_312Mbp.vcf.gz.
# (scripts for these steps not provided)


# Remove indels, keep multiallelic SNPs

bcftools view --exclude-types indels MorexV3_5H_67_312Mbp.vcf.gz > MorexV3_5H_67_312Mbp_noindels.vcf 
bgzip MorexV3_5H_67_312Mbp_noindels.vcf 

# Keep only domesticated barleys
selection="IPK_domesticated"

# Requires a txt file containing the names of accessions you want to keep. This can be easily done by filtering a passport data file by domestication status (not provided, can be found in Milner et al., 2019)
bcftools view -S ${selection}_names.txt MorexV3_5H_67_312Mbp_noindels.vcf.gz> ${selection}_5H.vcf 
bgzip ${selection}_5H.vcf

# Only keep SNPs with < 20% missing data and a minor allele count > 1
vcftools --gzvcf ${selection}_5H.vcf.gz --mac 1 --max-missing 0.8 --recode --recode-INFO-all --out ${selection}_5H_20missing_MAC 

# Convert the vcf subsets to hapmap
run_pipeline.pl -Xms10g -Xmx50g -fork1 -vcf ${selection}_5H_20missing_MAC.recode.vcf -export ${selection}_5H_20missing_MAC.hmp.txt -exportType Hapmap -runfork1


