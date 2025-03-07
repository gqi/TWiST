rm(list=ls())
library(dplyr)
library(plink2R)
library(TWiST)
sumstats <- data.table::fread("example_data/RA_sumstats_chr6.txt")
ngwas <- 22350*74823/(22350+74823)

# Load weights
ctype <- "T_CD8"
load(paste0("pretrained_models/twist_weights_",ctype,".rda"))
wgtlist.chr <- wgtlist %>% filter(CHR==6)
weights_pred.chr <- weights_pred[wgtlist.chr$ID]
genos.chr <- read_plink("example_data/1000G.EUR.6")

res <- twist_association(
    sumstat=sumstats, wgtlist=wgtlist.chr, weights_pred=weights_pred.chr,
    bim_train=bim_train, genos=genos.chr, ngwas=ngwas)

names(res)
str(res$out.tbl)

library(qqman)
png(filename=paste0("example_data/QQ_",ctype,"_chr6.png"), width=2500, height = 900, res=300)
par(mfrow=c(1,3))
qq(res$out.tbl$p.global, main="Global test", ylim=c(0,220))
qq(res$out.tbl$p.dynamic, main="Dynamic test", ylim=c(0,220))
qq(res$out.tbl$p.nonlinear, main="Nonlinear test", ylim=c(0,220))
dev.off()
