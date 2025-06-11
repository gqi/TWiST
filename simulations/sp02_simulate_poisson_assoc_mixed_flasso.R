# Association analysis using mixed models
rm(list=ls())
source("../../code/functions/sctwas_assoc_mixed.R")
ncausalsnpsvec <- c(5,10,20)
nknotssimvec <- c(1,3,5)
simmodvec <- c("null","unimodal","reverse","constant")
ngwas <- 1e5
degree <- 3

ptseq <- seq(0,1,by=0.01)
exon.all <- data.frame()

# Analysis model - set up B-spline basis
knotsana.w <- c(0.25, 0.5,0.75) # Number of B-spline knots for weights
knotsana.beta <- seq(0.05,0.95,by=0.05)
ptbsana.beta <- bs(ptseq, knots=knotsana.beta, degree=degree, Boundary.knots=c(0,1), intercept=TRUE)
omegalist.ana <- generate_omega(
    knots.w=knotsana.w, degree.w=degree, 
    knots.beta=knotsana.beta, degree.beta=degree)

# Obtain environment variables
temp <- as.integer(commandArgs(trailingOnly = TRUE)) # ncausalsnps (1:3), nknotssim (1:3), repid (1:5)
print(temp)

# Set simulation parameters based on environmental variables
ncausalsnpsval <- ncausalsnpsvec[temp[1]]
nknotssimval <- nknotssimvec[temp[2]]
repid <- temp[3] # 1:5
simmod <- simmodvec[temp[4]]
set.seed(35423*prod(temp))

# Load trained prediction model
filename <- paste0("../sp01_simulate_poisson_pred_flasso/res_ncausalsnps",ncausalsnpsval,"_nknotssim",nknotssimval,"_rep",repid,".rda")
load(filename)

ind <- sapply(weightlist, function(x) sum(x!=0)>0)
Wlist <- Wlist[ind]
exon <- exon[ind,]
weightlist <- weightlist[ind]
Rlist <- Rlist[ind]
weightlist.pb <- weightlist.pb[ind]

# Set up data frame to store association results
exon.temp <- exon %>% dplyr::select(gene=gene_id_sh, ncausalsnps, nknotssim) %>%
    mutate(p_global=NA, p_dynamic=NA, p_nonlinear=NA)
pbmethodvec <- c("lasso","ridge","enet")
for (pbmethod in pbmethodvec){
    exon.temp[[paste0("p_global_pb_",pbmethod)]] <- NA
    exon.temp[[paste0("p_dynamic_pb_",pbmethod)]] <- NA
}

# Mixed model association analysis
assoc.gene <- vector("list", length=nrow(exon.temp))
for (i in 1:length(Wlist)){ 
    # Simulate GWAS summary stats
    knotssim <- seq(0,1,length.out=exon$nknotssim[i]+2)[-c(1,exon$nknotssim[i]+2)]
    bspl <- create.bspline.basis(norder=degree+1, breaks=c(0,knotssim,1))
    coefmat <- diag(rep(1,bspl$nbasis))
    fd.bspl <- fd(coefmat, basisobj=bspl)
    omega.sim <- inprod.bspline(fd.bspl,fd.bspl)
    
    # Generate effect size
    beta <- rep(0, ncol(omega.sim))
    if (simmod=="unimodal"){
        midind <- ceiling(length(beta)/2)
        beta[midind] <- 0.2
        beta[c(midind-1,midind+1)] <- 0.1
    } else if (simmod=="reverse"){
        midind <- ceiling(length(beta)/2)
        beta[midind-1] <- 0.3
        beta[midind+1] <- -0.3
    } else if (simmod=="constant"){
        beta[1:length(beta)] <- 0.05
    }
    
    # Simulate summary statistics
    b <- Rlist[[i]]%*%Wlist[[i]]%*%omega.sim%*%beta
    bhat <- mvrnorm(n=1, mu=b, Sigma=Rlist[[i]]/ngwas)
    
    # Association analysis using mixed model
    assoc.gene[[i]] <- assoc_onegene(
        wgt.matrix=weightlist[[i]], bhat=bhat, R=Rlist[[i]], ngwas=ngwas, 
        omegal=omegalist.ana$omegal, omega=omegalist.ana$omega, omega2=omegalist.ana$omega2, 
        method="ml", logsigma2vec=seq(-20,10,by=0.2), test=TRUE)
    
    # Save association tests -  single-cell
    exon.temp$p_global[i] <- assoc.gene[[i]]$p.global
    exon.temp$p_dynamic[i] <- assoc.gene[[i]]$p.dynamic
    exon.temp$p_nonlinear[i] <- assoc.gene[[i]]$p.nonlinear
    
    # Association testing - pseudobulk
    for (pbmethod in pbmethodvec){
        # Global
        wgt.matrix.pb <- as.vector(weightlist.pb[[i]][[pbmethod]])
        if (sum(wgt.matrix.pb!=0)>0){
            zstat <- as.numeric(crossprod(wgt.matrix.pb, bhat*sqrt(ngwas))/
                                    sqrt(crossprod(wgt.matrix.pb,Rlist[[i]]%*%wgt.matrix.pb)))
            exon.temp[[paste0("p_global_pb_",pbmethod)]][i] <- 2*pnorm(-abs(zstat))
        }
        
        # Dynamic
        weights <- cbind(as.vector(weightlist.pb[[i]][[paste0(pbmethod,".early")]]), 
                         as.vector(weightlist.pb[[i]][[paste0(pbmethod,".late")]]))
        if (sum(colSums(weights!=0)==0)==0){ # No all zero columns
            betahat.pb <- solve(crossprod(weights,Rlist[[i]]%*%weights), crossprod(weights,bhat))
            var.betahat.pb <- solve(crossprod(weights,Rlist[[i]]%*%weights))/ngwas
            contr.vec <- c(1,-1)
            zstat <- crossprod(contr.vec,betahat.pb)/sqrt(crossprod(contr.vec,var.betahat.pb%*%contr.vec))
            exon.temp[[paste0("p_dynamic_pb_",pbmethod)]][i] <- 2*pnorm(-abs(zstat))
        }
     }
}

save(exon.temp, assoc.gene, file=paste0("simresults_ncausalsnps",ncausalsnpsval,"_nknotssim",nknotssimval,"_rep",repid,"_",simmod,".rda"))
