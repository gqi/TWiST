# Compute scTWAS weights
# Request multi-core job: salloc --mem=30G -t  60:00 -p gqi-32c256g --ntasks 1 --cpus-per-task 16 --mem-per-cpu 4G
# setwd("../../simulations/s05_simulate_pt_unknown")
rm(list=ls())
library(tidyverse)
library(plink2R)
library(splines)
ncausalsnpsval <- 10 # Number of causal SNPs
nknotssimval <- 5 # Number of B-spline knots for the simulated SNP effects
degree <- 3 # Degree of B-spline basis (cubic)

chrnum <- 6 # Chromosome number
# Load list of genes. Subset to ncausalsnps, nknotssim
load("../../Postdoc_research/single_cell/data/Cuomo_ASE/geno/exon_annotation_gencode_v19.rda")
exon <- exon %>% filter(feature=="gene" & seqname==chrnum) %>%
    mutate(gene_id_sh=gsub("\\.[[0-9]+$","",gene_id),
           tss=ifelse(strand=="+",start,end))
# Load genotype data
geno <- read_plink("TWiST/example_data/1000G.EUR.6")
colnames(geno$bim) <- c("CHR","SNP","cM","BP","A1","A2")
rownames(geno$bed) <- gsub("^0:","",rownames(geno$bed))

# Simulation parameters
sigma2.snp <- 0.01
sigma2.w <- 0.01
ncells <- 5000
pt <- runif(ncells, 0, 1) # Pseudotime
pt.eff <- 1
bpthr <- 5e5 # Range of cis-region

# Create individual-cell correspondence
donorind <- sample(nrow(geno$fam), ncells, replace=TRUE)

# Extract cis-SNPs of a gene, use the HLA-A as an example
gene.info <- exon %>% filter(gene_name=="HLA-A") # geneid_short is unique. gene.info should have 1 row.
bim.gene <- geno$bim %>% filter(BP<gene.info$tss+bpthr & BP>gene.info$tss-bpthr)
geno.gene <- geno$bed[,bim.gene$SNP,drop=F]
geno.gene <- scale(geno.gene)

##  Simulate gene expression
# Simulate SNP effect size on gene expression
knotssim <- seq(0,1,length.out=nknotssimval+2)[-c(1,nknotssimval+2)]
ptbs <- bs(pt, knots=knotssim, degree=3, Boundary.knots=c(0,1), intercept=TRUE)

W <- matrix(0, nrow=ncol(geno.gene), ncol=ncol(ptbs))
indcausal <- sample(nrow(W), ncausalsnpsval)
W[indcausal,] <- matrix(rep(rnorm(length(indcausal),mean=0,sd=sqrt(sigma2.snp)),ncol(ptbs)),nrow=length(indcausal)) +
    rnorm(length(indcausal)*ncol(W), mean=0, sd=sqrt(sigma2.w))

linpred <- geno.gene %*% W
linpred <- linpred[donorind,]
gene_exp_i <- rowSums(linpred*ptbs) + pt.eff*pt 
gene_exp_i <- rpois(length(gene_exp_i),lambda=exp(gene_exp_i))

## Load and simulate covariates
# Use plink to compute genotype PCs: ../../../software/plink1_9/plink --bfile 1000G.EUR.6 --pca --out 1000G_EUR_Chr6_PC
geno.pc <- data.table::fread("TWiST/example_data/1000G_EUR_Chr6_PC.eigenvec")
colnames(geno.pc) <- c("FID","IID",paste0("genPC_",1:20))
covar <- data.frame(age=sample(20:70,size=nrow(geno.pc),replace=TRUE),
                    sex=rbinom(nrow(geno.pc),size=1,prob=0.5), geno.pc%>%dplyr::select(genPC_1:genPC_10))
expression.pc <- matrix(rnorm(ncells*10,mean=0,sd=1), nrow=ncells,
                        dimnames=list(NULL,paste0("PC_",1:10)))
covar <- as.matrix(cbind(covar[donorind,], expression.pc))

# Dupliate genotype to cell level
geno.cell <- geno.gene[donorind,]

# Simulate library size
libsize=sample(10000:100000, size=nrow(covar), replace=TRUE)

# Save simulated data
save(gene_exp_i, geno.cell, pt, libsize, covar, file="TWiST/example_data/example_data_stage1training.rda")

# Training prediction model
library(TWiST)
model <- twist_train_model(y=gene_exp_i, geno_cell=geno.cell, pt=pt, knots=c(0.25,0.5,0.75), 
                           degree=3, lambda=NULL, nlambda=10, libsize=libsize, covar=covar)
