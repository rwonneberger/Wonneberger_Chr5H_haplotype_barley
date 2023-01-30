# Wonneberger_Chr5H_haplotype_barley
Scripts used in the paper "Major chromosome 5H haplotype switch structures the European two-rowed spring barley germplasm of the past 190 years" by Wonneberger et al., 2023

The analyses in this paper rely on a number of files:
- All analyses of variant data of the European two-rowed spring barley panel is based on the hapmap file published in Schreiber et al., 2023. The genotypes are encoded as A, C, G, T. For some of the analyses, the data needs to be encoded as AA, CC, GG, TT, and for other analyses a numeric file is needed. Conversion to diploid encoding and to numeric format are done with the `Prepare_snpmatrix_input_files.R` script
- A metadata file containing information about cultivars, year of release, country of origin, haplotype and cluster membership. This can be generated from Supplemental Data 1
- Analyses of the IPK gene bank collection require a SNP matrix of all accessions. [Milner et al., 2019](https://www.nature.com/articles/s41588-018-0266-x) published a SNP matrix with SNPs mapped to Morex V1. We used the same matrix but with SNPs mapped to Morex V3, provided by Dr. Martin Mascher, IPK. In [Milner et al., 2019](https://doi.ipk-gatersleben.de/DOI/ecfbdb3d-4882-406c-9e82-7758ed5395c7/4f58176f-4824-4c32-bca1-3d87500d82f3/2) you will also find passport information which you will need.
