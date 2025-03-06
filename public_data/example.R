rm(list=ls())
library(dplyr)
library(plink2R)
library(splines)
library(fda)
sumstats <- data.table::fread("public_data/RA_sumstats_chr6.txt")
ngwas <- 22350*74823/(22350+74823)

# Load weights
ctype <- "T_CD8"
load(paste0("public_data/twist_weights_",ctype,".rda"))
wgtlist$nsnps <- sapply(weights_pred, function(x) sum(rowSums(x$Wmat!=0)>0))
wgtlist$num_bsbases <- sapply(weights_pred, function(x) ncol(x$Wmat))
wgtlist.chr <- wgtlist %>% filter(CHR==6)
weights_pred.chr <- weights_pred[wgtlist.chr$ID]
genos.chr <- read_plink("public_data/1000G.EUR.6")

res <- twist_association(
    sumstat=sumstats, wgtlist=wgtlist.chr, weights_pred=weights_pred.chr,
    bim_train=bim_train, genos=genos.chr, ngwas=ngwas)
