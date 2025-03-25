# TWiST

TWiST (**TW**AS **i**n p**S**eudo**T**ime) is an R package for single-cell TWAS analysis of heterogeneous cell types, where gene expression and eQTL effects can vary along a continuous cell state within the cell type. 
Cell state is defined by pseudotime. This package implements two main analyses:

**Stage 1:** Train models to predict gene expression using cis-SNPs using single-cell eQTL data. Gene expression is modeled directly as count data. 

**Stage 2:**: Conduct association analysis between gene expression and trait. 

<img src="/example_data/overview.png" alt="overview" width="800"/>

**We have provided pre-trained models for CD4+ T cells, CD8+ T cells, and B cells. Users who are interested in these cell types may use the pre-trained models and skip Stage 1 (see example below).**

## Installation

Install `TWiST` from GitHUb:
```
devtools::install_github("gqi/TWiST")
```

In addition, install the `plink2R` package to read genotype data (PLINK files) into R:
```
devtools::install_github("gabraham/plink2R/plink2R")
```

## Pre-trained models

Download pre-trained models for three immune cell types from folder `pretrained_models`: `twist_weights_T_CD4.rda` (CD4+ T cells), `twist_weights_T_CD8.rda` (CD8+ T cells), `twist_weights_T_CD8.rda` (B cells).

Each `.rda` file includes three objects:
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
* `bim_train`: Information of model SNPs in `weights_pred` in the format of PLINK bim file. A data frame of the following columns:
    * `CHR`: Chromosome
    * `SNP`: SNP ID
    * `cM`: SNP position in centimorgan
    * `BP`: SNP position in base pair
    * `A1`: Effect allele. Coefficients in `weights_pred` are with respect to `A1`. <! -- PLINK bed file counts the number of A1 allele --> 
    * `A2`: Other allele

## Example: Association analysis
To run this example, download additional datasets provided in folder `example_data`. They include:
* `RA_sumstats_chr6.txt`: GWAS summary statistics for rheumatoid arthritis, chromosome 6 (Ishigaki et al, Nature Genetics 2022).
* `1000G.EUR.6.{bed,bim,fam}`: 1000 Genomes European genotype data, chromosome 6. Download all the chromosomes from [here](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2).

First, load required packages
```
library(dplyr)
library(plink2R)
library(TWiST)
```

Read GWAS summary statistics into R. Compute the effective sample size defined as ncases*ncontrols/(ncases+ncontrols).
```
sumstats <- data.table::fread("example_data/RA_sumstats_chr6.txt")
ngwas <- 22350*74823/(22350+74823)
```

Load pre-trained models and subset to chromosome 6 (using CD8+ T cells as an example) 
```
# Load weights
ctype <- "T_CD8"
load(paste0("pretrained_models/twist_weights_",ctype,".rda"))
wgtlist.chr <- wgtlist %>% filter(CHR==6)
weights_pred.chr <- weights_pred[wgtlist.chr$ID]
```

Load reference genotype data - 1000 Genomes European sample
```
genos.chr <- read_plink("example_data/1000G.EUR.6")
```

Run TWiST association analysis
```
res <- twist_association(
    sumstat=sumstats, wgtlist=wgtlist.chr, weights_pred=weights_pred.chr,
    bim_train=bim_train, genos=genos.chr, ngwas=ngwas)
```

View results (type `?twist_association` for definition of the outputs)
```
names(res)
# [1] "out.tbl"   "betal"     "var.betal" "beta"      "var.beta"  "knots"  

str(res$out.tbl)
# 'data.frame':	80 obs. of  10 variables:
#  $ ID         : chr  "ENSG00000112679" "ENSG00000170542" "ENSG00000124570" "ENSG00000214113" ...
#  $ CHR        : int  6 6 6 6 6 6 6 6 6 6 ...
#  $ P0         : int  291630 2887500 2948393 5102827 10723148 16129356 18224099 24705294 24804513 26104104 ...
#  $ P1         : int  351355 2903514 2972090 5261172 10731362 16148479 18265054 24721064 24936188 26104518 ...
#  $ tss        : int  291630 2903514 2972090 5261172 10723148 16129356 18265054 24721064 24936188 26104104 ...
#  $ sigma2     : num  2.06e-09 2.06e-09 2.06e-09 2.06e-09 2.06e-09 ...
#  $ p.global   : num  0.779 0.133 0.1534 0.7266 0.0599 ...
#  $ p.dynamic  : num  0.604 0.655 0.332 0.466 0.167 ...
#  $ p.nonlinear: num  0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 ...
#  $ degree     : num  3 3 3 3 3 3 3 3 3 3 ...
```

Create QQ plots for global, dynamic and nonlinear tests:
```
library(qqman)
par(mfrow=c(1,3))
qq(res$out.tbl$p.global, main="Global test", ylim=c(0,220))
qq(res$out.tbl$p.dynamic, main="Dynamic test", ylim=c(0,220))
qq(res$out.tbl$p.nonlinear, main="Nonlinear test", ylim=c(0,220))
```

<img src="/example_data/QQ_T_CD8_chr6.png" alt="QQ" width="800"/>

## Example: Training prediction models

If you have your own single-cell eQTL data and would like to train your own prediction model, below is an example using simulated data:

Load example dataset:
```
load("example_data_stage1training.rda")
```

This dataset includes real genotype data and genotype principal components (PCs) from the European subset of 1000 Genomes (chromosome 6) and simulated gene expression, librarize, and covariates:
* `gene_exp_i`: Expression count data for one gene.
* `geno_cell`: Genotype of cis-SNPs, duplicated to cell level (cells from the same individual have the same genotypes).
* `libsize`: Library size.
* `covar`: Covariates: age, sex, 10 expression PCs (`PC_{1:10}`), 10 genotype PCs (`genPC_{1:10}`).

Train prediction model (for one gene):
```
library(TWiST)
model <- twist_train_model(y=gene_exp_i, geno_cell=geno.cell, pt=pt, knots=c(0.25,0.5,0.75), 
                           degree=3, lambda=NULL, nlambda=10, libsize=libsize, covar=covar)
```

Aggregate models across genes into the format described in the previous section before proceeding to association analysis.

Codes for simulating this example dataset are provided [here](/example_data/simulate_example_training.R).

# Reference

Qi G, Lila E, Ji Z, Shojaie A, Battle A, Sun W. Transcriptome-wide association studies at cell state level using single-cell eQTL data. medRxiv (2025): 25324128. https://doi.org/10.1101/2025.03.17.25324128
