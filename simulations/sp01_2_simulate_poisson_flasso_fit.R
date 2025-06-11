# Compute scTWAS weights
rm(list=ls())
library(tidyverse)
library(plink2R)
library(splines)
library(fda)
library(glmnet)
library(parallel)
source("../../code/functions/sctwas_weights_poisson.R")
ncausalsnpsvec <- c(5,10,20)
nknotssimvec <- c(1,3,5)

temp <- as.integer(commandArgs(trailingOnly = TRUE)) # ncausalsnps (1:3), nknotssim (1:3), repid (1:3); batch (1:10)
print(temp)

ncausalsnpsval <- ncausalsnpsvec[temp[1]]
nknotssimval <- nknotssimvec[temp[2]]
repid <- temp[3]
batch <- temp[4]
set.seed(34254*prod(temp))

load(paste0("geneexp_donorind_ncausalsnps",ncausalsnpsval,"_nknotssim",nknotssimval,"_rep",repid,".rda"))
batchind <- (20*(batch-1)+1):(20*batch)
exon <- exon[batchind,]
res <- res[batchind] 
geno.cis.split <- geno.cis.split[batchind] 

# Initialize storage lists
Wlist <- Rlist <- weightlist <- lambdavec <- 
    weightlist.pb <- vector("list",length=length(res))
outname <- paste0("res_ncausalsnps",ncausalsnpsval,"_nknotssim",nknotssimval,"_rep",repid,"_batch",batch,".rda")
for (i in 1:length(geno.cis.split)){
    print(i)
    geno.dat <- geno.cis.split[[i]]
    geno.cell <- geno.dat$geno.gene[donorind,]
    y <- geno.dat$gene_exp_i
    
    desmat <- generate_desmat(geno.cell, ptbsana)
    grp <- rep(1:ncol(geno.cell),each=ncol(ptbsana))
    
    # Cross validation to select optimal lambda parameter
    gr_cv <- cv.grpreg(X=desmat, y=y, group=grp, penalty="grLasso",
                       family="poisson", nlambda=20, alpha=0.5, nfolds=5)
    res.opt <- grpreg(X=desmat, y=y, group=grp, penalty="grLasso",
                      family="poisson", lambda=gr_cv$lambda.min, alpha=0.5)
    # Different from glmnet, res.opt$beta from grpnet includes intercept
    betamat.lm <- matrix(res.opt$beta[-1], nrow=ncol(geno.cell), byrow=TRUE,
                         dimnames=list(colnames(geno.cell),paste0("bs",1:ncol(ptbsana))))
    
    # Pseudobulk prediction models
    # Global test
    df <- data.frame(donor=donorind, y_sc=y) %>% group_by(donor) %>%
        summarize(y=log(1+mean(y_sc))) %>%
        arrange(donor)
    
    cv <- cv.glmnet(x=geno.dat$geno.gene[df$donor,], y=df$y, family="gaussian", alpha=1)
    mod.lasso <- glmnet(x=geno.dat$geno.gene[df$donor,], y=df$y, family="gaussian", alpha=1, lambda=cv$lambda.min)
    cv <- cv.glmnet(x=geno.dat$geno.gene[df$donor,], y=df$y, family="gaussian", alpha=0)
    mod.ridge <- glmnet(x=geno.dat$geno.gene[df$donor,], y=df$y, family="gaussian", alpha=0, lambda=cv$lambda.min)
    cv <- cv.glmnet(x=geno.dat$geno.gene[df$donor,], y=df$y, family="gaussian", alpha=0.5)
    mod.enet <- glmnet(x=geno.dat$geno.gene[df$donor,], y=df$y, family="gaussian", alpha=0.5, lambda=cv$lambda.min)
    
    # Dynamic test
    dfearly <- data.frame(donor=donorind, y_sc=y, ptest) %>% 
        filter(ptest<0.5) %>%
        group_by(donor) %>%
        summarize(y=log(1+mean(y_sc))) %>%
        arrange(donor)
    dflate <- data.frame(donor=donorind, y_sc=y, ptest) %>% 
        filter(ptest>=0.5) %>%
        group_by(donor) %>%
        summarize(y=log(1+mean(y_sc))) %>%
        arrange(donor)
    
    cv <- cv.glmnet(x=geno.dat$geno.gene[dfearly$donor,], y=dfearly$y, family="gaussian", alpha=1)
    mod.lasso.early <- glmnet(x=geno.dat$geno.gene[dfearly$donor,], y=dfearly$y, family="gaussian", alpha=1, lambda=cv$lambda.min)
    cv <- cv.glmnet(x=geno.dat$geno.gene[dfearly$donor,], y=dfearly$y, family="gaussian", alpha=0)
    mod.ridge.early <- glmnet(x=geno.dat$geno.gene[dfearly$donor,], y=dfearly$y, family="gaussian", alpha=0, lambda=cv$lambda.min)
    cv <- cv.glmnet(x=geno.dat$geno.gene[dfearly$donor,], y=dfearly$y, family="gaussian", alpha=0.5)
    mod.enet.early <- glmnet(x=geno.dat$geno.gene[dfearly$donor,], y=dfearly$y, family="gaussian", alpha=0.5, lambda=cv$lambda.min)
    
    cv <- cv.glmnet(x=geno.dat$geno.gene[dflate$donor,], y=dflate$y, family="gaussian", alpha=1)
    mod.lasso.late <- glmnet(x=geno.dat$geno.gene[dflate$donor,], y=dflate$y, family="gaussian", alpha=1, lambda=cv$lambda.min)
    cv <- cv.glmnet(x=geno.dat$geno.gene[dflate$donor,], y=dflate$y, family="gaussian", alpha=0)
    mod.ridge.late <- glmnet(x=geno.dat$geno.gene[dflate$donor,], y=dflate$y, family="gaussian", alpha=0, lambda=cv$lambda.min)
    cv <- cv.glmnet(x=geno.dat$geno.gene[dflate$donor,], y=dflate$y, family="gaussian", alpha=0.5)
    mod.enet.late <- glmnet(x=geno.dat$geno.gene[dflate$donor,], y=dflate$y, family="gaussian", alpha=0.5, lambda=cv$lambda.min)
    
    weight.pb <- list(lasso=mod.lasso$beta, ridge=mod.ridge$beta, enet=mod.enet$beta,
                      lasso.early=mod.lasso.early$beta, ridge.early=mod.ridge.early$beta, enet.early=mod.enet.early$beta,
                      lasso.late=mod.lasso.late$beta, ridge.late=mod.ridge.late$beta, enet.late=mod.enet.late$beta)
    
    Wlist[[i]] <- res[[i]]$W
    Rlist[[i]] <- res[[i]]$R
    weightlist[[i]] <- betamat.lm
    lambdavec[[i]] <- gr_cv$lambda.min
    weightlist.pb[[i]] <- weight.pb
    
    save(Wlist, Rlist, omegalist, omega2list, weightlist, weightlist.pb, lambdavec, exon, ptmat, file=outname)
}

