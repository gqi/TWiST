# Assemble prediction models
rm(list=ls())
library(tidyverse)
ncausalsnpsvec <- c(5,10,20)
nknotssimvec <- c(1,3,5)

for (ncausalsnpsval in ncausalsnpsvec){
    for (nknotssimval in nknotssimvec){
        for (repid in 1:3){
            print(c(ncausalsnpsval, nknotssimval, repid))
            
            Wlist.all <- Rlist.all <- weightlist.all <- weightlist.pb.all <- lambdavec.all <- list()
            exon.all <- data.frame()
            for (batch in 1:10){
                outname <- paste0("res_ncausalsnps",ncausalsnpsval,"_nknotssim",nknotssimval,"_rep",repid,"_batch",batch,".rda")
                load(outname)
                
                Wlist.all <- c(Wlist.all, Wlist)
                Rlist.all <- c(Rlist.all, Rlist)
                weightlist.all <- c(weightlist.all, weightlist)
                weightlist.pb.all <- c(weightlist.pb.all, weightlist.pb)
                lambdavec.all <- c(lambdavec.all, lambdavec)
                exon.all <- rbind(exon.all, exon)
            }
            
            Wlist <- Wlist.all
            Rlist <- Rlist.all
            weightlist <- weightlist.all
            weightlist.pb <- weightlist.pb.all
            lambdavec <- lambdavec.all
            exon <- exon.all
            
            filename <- paste0("res_ncausalsnps",ncausalsnpsval,"_nknotssim",nknotssimval,"_rep",repid,".rda")
            save(Wlist, Rlist, omegalist, omega2list, weightlist, weightlist.pb, lambdavec, exon, ptmat, file=filename)
            rm(Wlist, Rlist, omegalist, omega2list, weightlist, weightlist.pb, lambdavec, exon, ptmat, filename, outname,
               Wlist.all, Rlist.all, weightlist.all, weightlist.pb.all, lambdavec.all, exon.all)
        }
    }
}


