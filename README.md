# TWiST

TWiST (**TW**AS **i**n p**S**eudo**T**ime) is an R package for single-cell TWAS analysis of heterogeneous cell types, where gene expression and eQTL effects can vary along a continuous cell state within the cell type. 
Cell state is defined by pseudotime. 

Stage 1 of TWiST is to train a model to predict single-cell gene expression using single nucleotide polymorphisms (SNPs) in the cis-region of the gene. 
TWiST addresses two main challenges unique to single-cell data: 1) The data is highly sparse and cannot be modeled using normal distribution (as in bulk TWAS) even after transformation. 2) There are many cells per individual which all have the same genotype, hence standard sparse regression methods cannot be applied. 

To resolve the challenges, TWiST models gene expression counts directly using Poisson distribution. 
The log-transformed mean of the Poisson is linked to a continuous function along pseudotime which represents the genetically regulated gene expression (GReX).

In Stage 2, we model the trait (y_i) by the aggregated effect of GReX across pseudotime t: y_i=d_0+∫β(t) v_i (t)dt+e_i. The effect of GReX on trait β(t) is modeled as a spline function along pseudotime. This model can be fitted using individual-level data or GWAS summary statistics. Coefficient β(t) can be used for downstream hypothesis testing of global (existence of any effect at any pseudotime point), dynamic (effect varies over pseudotime), and nonlinear effects. Details are described in Methods. 

We have provided pre-trained models for CD4+ T cells, CD8+ T cells, and B cells. Users who are interested in conducting TWAS for these cell types may skip Step 1. Step-by-step tutorial is shown below

### Installation
```
devtools::install_github("gqi/TWiST")
```

### Example

To run this example, download the data provided in folder `public_data`. They include
* `twist_weights_{T_CD4,T_CD8,B}.rda`: Pre-trained model for CD4+ T cells, CD8+ T cells, and B cells in the OneK1K data.
* `RA_sumstats_chr6.txt`: GWAS summary statistics for rheumatoid arthritis, chromosome 6 (Ishigaki et al, Nature Genetics 2022).
* `1000G.EUR.6.{bed,bim,fam}`: 1000 Genomes European genotype data, chromosome 6. Download all the chromosomes from [here](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2).
