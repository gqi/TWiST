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

Install `TWiST` from GitHUb:
```
devtools::install_github("gqi/TWiST")
```

In addition, install the `plink2R` package to read genotype data (PLINK files) into R:
```
devtools::install_github("gabraham/plink2R/plink2R")
```

### Example

To run this example, download the data provided in folder `public_data`. They include
* `twist_weights_{T_CD4,T_CD8,B}.rda`: Pre-trained model for CD4+ T cells, CD8+ T cells, and B cells in the OneK1K data.
* `RA_sumstats_chr6.txt`: GWAS summary statistics for rheumatoid arthritis, chromosome 6 (Ishigaki et al, Nature Genetics 2022).
* `1000G.EUR.6.{bed,bim,fam}`: 1000 Genomes European genotype data, chromosome 6. Download all the chromosomes from [here](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2).

Load pre-trained models, using CD8+ T cells as an example:
```
ctype <- "T_CD8"
load(paste0("public_data/twist_weights_",ctype,".rda"))
```

This `.rda` file includes three objects:
* `wgtlist`: Information of genes for which the model has been trained. A data frame of five columns:
    * `ID`: Gene ID
    * `CHR`: Chromosome
    * `P0`: Gene start
    * `P1`: Gene end
    * `tss`: Transcription start site
* `weights_pred`: Pre-trained prediction models. A list of which each entry is a gene, in the same order as `wgtlist$ID`. For each genes, the following four entries are available:
    * `Wmat`: Coefficients of pre-trained prediction models. A matrix of (number of SNPs) x (number of B-spline bases).
    * `knots`: Internal knots for B-spline basis functions that are used to model model SNP effects on gene expression. (boundary knots 0 and 1 are not included).
    * `degree`: Degree of B-spline basis functions. 
    * `n`: Number of cells in the eQTL data for model training.
* `bim_train`: PLINK bim file from the training data where prediction models were built. Should include all the SNPs in `weights_pred`. A data frame including the following columns
    * `CHR`: Chromosome
    * `SNP`: SNP ID
    * `cM`: SNP position in centimorgan
    * `BP`: SNP position in base pair
    * `A1`: Effect allele
    * `A2`: Other allele



