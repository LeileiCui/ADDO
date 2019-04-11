# ADDO

Author: Leilei Cui and Bin Yang

A comprehensive toolkit to detect, classify and visualise additive and non-additive Quantitative Trait Loci

## Motivation
Additivity and dominance are the major genetic components underlying variations in complex traits. During the past decade, genome-wide association studies (GWAS) have been used to map quantitative trait loci (QTLs) underlying complex traits. However, most GWAS focus on additive genetic effects while ignoring non-additive effects, on the assumption that most QTL act additively. Consequently, QTLs driven by dominance and other non-additive effects could be overlooked.

## Results
We developed ADDO, a highly-efficient tool designed to detect, classify and visualize quantitative trait loci (QTLs) with additive and non-additive effects. ADDO implements a mixed-model transformation to control for population structure and unequal relatedness that accounts for both additive and dominant genetic covariance among individuals, and decomposes single nucleotide polymorphism (SNP) effects into additive, partial dominance, dominance and overdominance categories. A matrix multiplication approach is used to accelerate the computation: a genome scan on 20 million markers from 836 individuals takes about 8.5 hours with 10 CPUs.

## Prerequisites
The following command line tools:
* gcta v.1.26 https://cnsgenomics.com/software/gcta/gcta_1.92.1beta5.zip
* plink v.1.90 http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip

The following R packages:
* data.table
* parallel
* bigmemory
* mvtnorm (Only required by the Heterotic Model)
* MASS (Only required by the Heterotic Model)
* GenABEL (optional)
* emma (optional)

(Note emma and gcta are used to calcuate the kinship matrix so only one is required)

## Running the examples
The following example dataset are available in the `data` directory:
TEST.bed : The genotypes information
TEST.bim : The loci information
TEST.fam : The individuals information
TEST.phe : The phenotypes data
TEST.covs : The covariates data

Dominant effect detection:
```
R CMD BATCH demo/TEST_AddDom_Model.r
```

Over-dominant effect detection:
```
R CMD BATCH demo/TEST_Heterotic_Model.r
```

## Pipeline


## Contact
leileicui_xuan@hotmail.com
