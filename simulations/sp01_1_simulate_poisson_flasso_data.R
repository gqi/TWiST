# Simulate gene expression data
rm(list=ls())
library(tidyverse)
library(plink2R)
library(splines)
library(fda)
library(glmnet)
source("../../code/functions/sctwas_weights_poisson.R")
ncausalsnpsvec <- c(5,10,20)
nknotssimvec <- c(1,3,5)

temp <- as.integer(commandArgs(trailingOnly = TRUE)) # ncausalsnps (1:3), nknotssim (1:3), repid (1:3)
ncausalsnpsval <- ncausalsnpsvec[temp[1]]
nknotssimval <- nknotssimvec[temp[2]]
repid <- temp[3]
set.seed(34254*prod(temp))

# Load list of genes. Subset to ncausalsnps, nknotssim
load("../simulate_expr_trait_onekgeno/exon_sim.rda")
exon <- exon %>% filter(ncausalsnps==ncausalsnpsval & nknotssim==nknotssimval)

# Load genotype data
geno <- read_plink("../../data/onek1k/geno_imputed_hm3/chr1_maf01hm3")
colnames(geno$bim) <- c("CHR","SNP","cM","BP","A1","A2")
rownames(geno$bed) <- gsub("^0:","",rownames(geno$bed))

# Analysis parameters
bpthr <- 5e5
degree <- 3
ngwas <- 1e5
nknotsana <- 3

# Simulation parameters
sigma2.snp <- 0.01
sigma2.w <- 0.01
ncells <- 50000
pt <- runif(ncells, 0, 1) # Pseudotime
pt.eff <- 1

# Compute proportion of expression variance explained by SNPs
exon$h2exp <- sapply(1:nrow(exon), function(i) {
    knotssim <- seq(0,1,length.out=exon$nknotssim[i]+2)[-c(1,exon$nknotssim[i]+2)]
    ptbs <- bs(pt, knots=knotssim, degree=3, Boundary.knots=c(0,1), intercept=TRUE)
    return(exon$ncausalsnps[i]*(sigma2.snp+sigma2.w)*mean(rowSums(ptbs^2)))
})

gene_exp <- matrix(0, nrow=nrow(exon), ncol=ncells)
Rlist <- Wlist <- list()
omegalist <- omega2list <- list()

knotsana <- seq(0,1,length.out=nknotsana+2)[-c(1,nknotsana+2)]
ptbsana <- bs(pt, knots=knotsana, degree=3, Boundary.knots=c(0,1), intercept=TRUE)
bspl <- create.bspline.basis(norder=degree+1, breaks=c(0,knotsana,1))
coefmat <- diag(rep(1,bspl$nbasis))
fd.bspl <- fd(coefmat, basisobj=bspl)

omegalist[[paste0("nknotsana_",nknotsana)]] <- inprod.bspline(fd.bspl,fd.bspl)
omega2list[[paste0("nknotsana_",nknotsana)]] <- inprod.bspline(fd.bspl,fd.bspl, nderiv1=2, nderiv2=2)

# Create individual-cell correspondence
donorind <- sample(nrow(geno$fam), ncells, replace=TRUE)

# Split genotype data into chunks for parallel computing
system.time(geno.cis.split <- lapply(1:nrow(exon), function(i) {
    gene <- exon$gene_id_sh[i]
    gene.info <- exon[i,] # geneid_short is unique. gene.info should have 1 row.
    bim.gene <- geno$bim %>% filter(BP<gene.info$tss+bpthr & BP>gene.info$tss-bpthr)
    geno.gene <- geno$bed[,bim.gene$SNP,drop=F]
    geno.gene <- scale(geno.gene)
    
    return(list(geno.gene=geno.gene, # bim.gene=bim.gene, 
                nknotssim=exon$nknotssim[i], ncausalsnps=exon$ncausalsnps[i]))
}))
names(geno.cis.split) <- exon$gene_id_sh

# Need to request cores for the interactive jobs. Otherwise one core will be allocated, no acceleration will be achieved.
system.time(res <- lapply(geno.cis.split, function(geno.dat){ # nrow(exon)
    geno.gene <- geno.dat$geno.gene
    
    knotssim <- seq(0,1,length.out=geno.dat$nknotssim+2)[-c(1,geno.dat$nknotssim+2)]
    ptbs <- bs(pt, knots=knotssim, degree=3, Boundary.knots=c(0,1), intercept=TRUE)
    
    W <- matrix(0, nrow=ncol(geno.gene), ncol=ncol(ptbs))
    indcausal <- sample(nrow(W), geno.dat$ncausalsnps)
    
    W[indcausal,] <- 
        matrix(rep(rnorm(length(indcausal),mean=0,sd=sqrt(sigma2.snp)),ncol(ptbs)),nrow=length(indcausal)) +
        rnorm(length(indcausal)*ncol(W), mean=0, sd=sqrt(sigma2.w))
    
    linpred <- geno.gene %*% W
    linpred <- linpred[donorind,]
    gene_exp_i <- rowSums(linpred*ptbs) + pt.eff*pt # + rnorm(ncells,0,sqrt(1-exon$h2exp[i]))
    
    # Simulate GWAS summary stats
    R <- cor(geno.gene)
    
    return(list(W=W, R=R, gene_exp_i=gene_exp_i))
}))

# Simulate expression using poisson
gene_exp <- sapply(res, function(x) rpois(length(x$gene_exp_i),lambda=exp(x$gene_exp_i)))
pcs <- log(1+gene_exp) %*% prcomp(log(1+gene_exp))$rotation

# Learn pseudotime
if (cor(pcs[,1], pt)>0){
    ptest <- pcs[,1]
} else{
    ptest <- -pcs[,1]
}
ptest <- rank(ptest)/length(ptest)
ptmat <- cbind(ptest, pt)
head(ptmat)
print(cor(pt,ptest))

# Assemble dataset
for (i in 1:length(geno.cis.split)){
    geno.cis.split[[i]]$gene_exp_i <- gene_exp[,i]
}

knotsana <- seq(0,1,length.out=nknotsana+2)[-c(1,nknotsana+2)]
ptbsana <- bs(ptest, knots=knotsana, degree=3, Boundary.knots=c(0,1), intercept=TRUE)
save(exon, ptmat, ptest, res, geno.cis.split, gene_exp, donorind, ptbsana, knotsana, omegalist, omega2list,
     file=paste0("geneexp_donorind_ncausalsnps",ncausalsnpsval,"_nknotssim",nknotssimval,"_rep",repid,".rda"))
