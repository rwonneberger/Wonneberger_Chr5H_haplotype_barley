#!/bin/bash

########################################################################
## Calculate FST between k-means clusters in windows along the genome ##
########################################################################

# This script will require the following programs: Tassel, vcftools, tabix.
# For the publication Tassel 5.2.49 and vcftools 0.1.16 were used

vcf_name="Hapmap_haploid"

##### Prepare datasets
run_pipeline.pl -Xmx5g -fork1 -h ${vcf_name}.hmp.txt -export -exportType VCF -runfork1

gzip -c ${vcf_name}.vcf > ${vcf_name}.vcf.gz

##### Parameters
window_name="100kb"
step_name="10kb"
size=100000
step=10000


##### Run FST analysis between clusters

vcftools --gzvcf ${vcf_name}.vcf.gz --weir-fst-pop Pop_DAPC_k2_Cluster1.txt --weir-fst-pop Pop_DAPC_k2_Cluster2.txt --fst-window-size ${size} --fst-window-step ${step} --out FST_DAPC_k2_C1_C2_${window_name}_window_${step_name}_step

vcftools --gzvcf ${vcf_name}.vcf.gz --weir-fst-pop Pop_DAPC_k3_Cluster1.txt --weir-fst-pop Pop_DAPC_k3_Cluster2.txt --fst-window-size ${size} --fst-window-step ${step} --out FST_DAPC_k3_C1_C2_${window_name}_window_${step_name}_step
vcftools --gzvcf ${vcf_name}.vcf.gz --weir-fst-pop Pop_DAPC_k3_Cluster1.txt --weir-fst-pop Pop_DAPC_k3_Cluster3.txt --fst-window-size ${size} --fst-window-step ${step} --out FST_DAPC_k3_C1_C3_${window_name}_window_${step_name}_step
vcftools --gzvcf ${vcf_name}.vcf.gz --weir-fst-pop Pop_DAPC_k3_Cluster2.txt --weir-fst-pop Pop_DAPC_k3_Cluster3.txt --fst-window-size ${size} --fst-window-step ${step} --out FST_DAPC_k3_C2_C3_${window_name}_window_${step_name}_step

