## QC alleles from summary stats: flipping strand, etc
allele.qc = function(a1,a2,ref1,ref2) {
    a1 = toupper(a1)
    a2 = toupper(a2)
    ref1 = toupper(ref1)
    ref2 = toupper(ref2)

    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip

    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;

    snp = list()
    snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
    snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
    snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

    return(snp)
}

generate_omega <- function(knots.w, degree.w, knots.beta, degree.beta){
    pteval <- seq(0.0005,0.9995,by=0.001)
    ptbs.w <- bs(pteval, knots=knots.w, degree=degree.w, Boundary.knots=c(0,1), intercept=TRUE)
    omegal <- 0.001*crossprod(ptbs.w,cbind(1,pteval))

    # Set up B-splines bases for W
    bspl.w <- create.bspline.basis(norder=degree.w+1, breaks=c(0,knots.w,1))
    fd.bspl.w <- fd(coef=diag(rep(1,bspl.w$nbasis)), basisobj=bspl.w)
    # B-spline bases for beta
    bspl.beta <- create.bspline.basis(norder=degree.beta+1, breaks=c(0,knots.beta,1))
    fd.bspl.beta <- fd(coef=diag(rep(1,bspl.beta$nbasis)), basisobj=bspl.beta)
    # Compute inner product
    omega <- inprod(fd.bspl.w,fd.bspl.beta)
    omega <- omega[,-c(1,ncol(omega))]
    omega2 <- inprod(fd.bspl.beta,fd.bspl.beta, Lfdobj1=2, Lfdobj2=2)
    omega2 <- omega2[-c(1,nrow(omega2)),-c(1,ncol(omega2))]

    return(list(omegal=omegal, omega=omega, omega2=omega2))
}

# QC for each gene
qc_impute_onegene <- function(wgt, weights_pred, wgtlist, w, bim_train, genos, sumstat, opt){

    wgt.matrix <- wgt[[w]]

    # Compute inner product of B-spline bases
    omegalist <- generate_omega(
        knots.w=weights_pred[[w]]$knots, degree.w=weights_pred[[w]]$degree,
        knots.beta=opt$knots_beta, degree.beta=opt$degree_beta)

    # Remove NAs (these should not be here)
    wgt.matrix[is.na(wgt.matrix)] = 0

    # Match up the SNPs and weights
    snps = bim_train[match(rownames(wgt.matrix),bim_train$SNP),] # New
    class(snps) <- "data.frame"
    m = match( snps[,2] , genos$bim[,2] )
    # Based on current setup, prediction models are trained among 1000G SNPs
    # Hence m.keep should be all TRUE
    m.keep = !is.na(m)
    snps = snps[m.keep,]
    wgt.matrix = wgt.matrix[m.keep,,drop=F]
    cur.genos = scale(genos$bed[,m[m.keep]])
    cur.bim = genos$bim[m[m.keep],]
    # Flip WEIGHTS for mismatching alleles between the weight matrix and the reference genome
    qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
    wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

    cur.FAIL = FALSE

    # Match up the SNPs and the summary stats
    m = match(cur.bim[,2] , sumstat$SNP)
    cur.Z = sumstat$Z[m]

    # Compute LD matrix
    cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)
    cur.miss = is.na(cur.Z)
    # Impute missing Z-scores
    if ( sum(cur.miss) != 0 ) { # If there are some missing SNPs
        if ( sum(!cur.miss) == 0 ) {
            cat( "WARNING : " , unlist(wgtlist[w,]) , "had no overlapping GWAS Z-scores, skipping this gene\n")
            cur.FAIL = TRUE
        } else if ( mean(cur.miss) > opt$max_impute ) {
            cat( "WARNING : " , unlist(wgtlist[w,]) , "had" , sum(cur.miss) , "/" , length(cur.miss) , "non-overlapping GWAS Z-scores, skipping this gene.\n")
            cur.FAIL = TRUE
        } else {
            cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
            cur.impz = cur.wgt %*% cur.Z[!cur.miss]
            cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
            cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)

            all.r2pred = rep(1,length(cur.Z))
            all.r2pred[ cur.miss ] = cur.r2pred
            if ( sum(is.na(all.r2pred)) != 0 ) {
                cat( "WARNING : " , unlist(wgtlist[w,]) , "had missing GWAS Z-scores that could not be imputed, skipping this gene.\n" )
                cur.FAIL = TRUE
            } else if ( mean( all.r2pred[ rowSums(wgt.matrix!=0)>0 ] ) < opt$min_r2pred ) {
                cat( "WARNING : " , unlist(wgtlist[w,]) , "had mean GWAS Z-score imputation r2 of" , mean( all.r2pred[rowSums(wgt.matrix!=0)>0] ) , "at expression weight SNPs, skipping this gene.\n")
                cur.FAIL = TRUE
            }
        }
    }

    if ( !cur.FAIL ) {
        return.list <- c(list(wgt.matrix=wgt.matrix, z=cur.Z, R=cur.LD), omegalist)
    } else{
        return(NULL)
    }
}

# log-likelihood of the equivalent mixed model
llkl_mixed <- function(logsigma2, betal, wgt.matrix, bhat, eigR, ngwas, omegal, omega, omega2){
    sigma2 <- exp(logsigma2)

    X.lin <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omegal
    X.nl <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omega
    bhat.decor <- crossprod(eigR$vectors,bhat)

    V <- diag(eigR$values)/ngwas + X.nl %*% solve(omega2, t(X.nl))*sigma2

    llkl <- -0.5*determinant(V,logarithm=TRUE)$modulus -
        0.5*crossprod(bhat.decor-X.lin%*%betal, solve(V,bhat.decor-X.lin%*%betal))

    return(as.numeric(llkl))
}

# Maximize the loglikelihood - all parameters
llkl_null_global <- function(wgt.matrix, bhat, eigR, ngwas, omegal, omega, omega2){

    llkl_mixed(logsigma2=-Inf, betal=c(0,0), wgt.matrix=wgt.matrix, bhat=bhat,
               eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2)
}

llkl_null_dynamic <- function(wgt.matrix, bhat, eigR, ngwas, omegal, omega, omega2){

    X.lin <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omegal[,1]
    bhat.decor <- crossprod(eigR$vectors,bhat)
    V <- diag(eigR$values)/ngwas
    betal0 <- solve(t(X.lin)%*%solve(V,X.lin),t(X.lin)%*%solve(V,bhat.decor))

    llkl_mixed(logsigma2=-Inf, betal=c(betal0,0), wgt.matrix=wgt.matrix, bhat=bhat,
               eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2)
}

llkl_null_nonlinear <- function(wgt.matrix, bhat, eigR, ngwas, omegal, omega, omega2){

    X.lin <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omegal
    bhat.decor <- crossprod(eigR$vectors,bhat)
    V <- diag(eigR$values)/ngwas
    betal <- solve(t(X.lin)%*%solve(V,X.lin),t(X.lin)%*%solve(V,bhat.decor))

    llkl_mixed(logsigma2=-Inf, betal=betal, wgt.matrix=wgt.matrix, bhat=bhat,
               eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2)
}

# REML criterion function
reml_mixed <- function(logsigma2, wgt.matrix, bhat, eigR, ngwas, omegal, omega, omega2){
    sigma2 <- exp(logsigma2)

    X.lin <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omegal
    X.nl <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omega
    bhat.decor <- crossprod(eigR$vectors,bhat)

    V <- diag(eigR$values)/ngwas + X.nl %*% solve(omega2, t(X.nl))*sigma2
    betal <- solve(t(X.lin)%*%solve(V,X.lin),t(X.lin)%*%solve(V,bhat.decor))

    reml <- -0.5*determinant(V,logarithm=TRUE)$modulus -
        0.5*crossprod(bhat.decor-X.lin%*%betal, solve(V,bhat.decor-X.lin%*%betal)) -
        0.5*determinant(crossprod(X.lin,solve(V,X.lin)),logarithm=TRUE)$modulus

    return(as.numeric(reml))
}

ml_mixed <- function(logsigma2, wgt.matrix, bhat, eigR, ngwas, omegal, omega, omega2){
    sigma2 <- exp(logsigma2)

    X.lin <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omegal
    X.nl <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omega
    bhat.decor <- crossprod(eigR$vectors,bhat)

    V <- diag(eigR$values)/ngwas + X.nl %*% solve(omega2, t(X.nl))*sigma2
    betal <- solve(t(X.lin)%*%solve(V,X.lin),t(X.lin)%*%solve(V,bhat.decor))

    ml <- -0.5*determinant(V,logarithm=TRUE)$modulus -
        0.5*crossprod(bhat.decor-X.lin%*%betal, solve(V,bhat.decor-X.lin%*%betal))

    return(as.numeric(ml))
}

## Likelihood ratio test based on random-effects model
# Here omega l is the inner product matrix for linear terms
# omega and omega2 do not include linear terms
assoc_onegene <- function(
        wgt.matrix, bhat, R, ngwas, omegal, omega, omega2, method="ml",
        logsigma2vec=seq(-20,10,by=0.2), test=FALSE){

    # Check whether R is symmetric, if not, make it symmetric
    if (!isSymmetric(R)){
        R <- (R+t(R))/2
    }

    rk <- qr(R)$rank

    if (rk>=5){

        eigR <- eigen(R)
        eigR$values <- eigR$values[1:rk]
        eigR$vectors <- eigR$vectors[,1:rk, drop=FALSE]
        # Transform predictors and summary stats
        X.lin <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omegal
        X.nl <- diag(eigR$values)%*%t(eigR$vectors)%*%wgt.matrix%*%omega
        bhat.decor <- crossprod(eigR$vectors,bhat)

        # Check conditions
        stopifnot(ncol(wgt.matrix)==nrow(omega) & nrow(wgt.matrix)==nrow(R))

        # Estimate variance component using REML
        if (method=="reml"){
            remlvec <- sapply(
                logsigma2vec, reml_mixed, wgt.matrix=wgt.matrix, bhat=bhat,
                eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2)
            sigma2 <- exp(logsigma2vec[which.max(remlvec)])

        } else if (method=="ml"){
            mlvec <- sapply(
                logsigma2vec, ml_mixed, wgt.matrix=wgt.matrix, bhat=bhat,
                eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2)
            sigma2 <- exp(logsigma2vec[which.max(mlvec)])
        }

        # Update betal with closed form
        V <- diag(eigR$values)/ngwas + X.nl %*% solve(omega2, t(X.nl))*sigma2
        betal <- solve(t(X.lin)%*%solve(V,X.lin),t(X.lin)%*%solve(V,bhat.decor))
        var.betal <- solve(crossprod(X.lin, solve(V,X.lin)))

        # BLUP
        C <- sigma2*solve(omega2, t(X.nl))
        betahat <- as.vector(C%*%solve(V,bhat.decor-X.lin%*%betal))
        # Variance of BLUP
        hatmat <- diag(1,nrow=length(bhat.decor)) -
            X.lin%*%solve(t(X.lin)%*%solve(V,X.lin),t(X.lin)%*%solve(V))
        var.betahat <- C%*%solve(V,hatmat)%*%V%*%t(hatmat)%*%solve(V,t(C))

        return.list <- list(betal=as.vector(betal), var.betal=var.betal, sigma2=sigma2,
                            betahat=betahat, var.betahat=var.betahat) #

        if (test){
            llkl.full <- llkl_mixed(
                logsigma2=log(sigma2), betal=betal, wgt.matrix=wgt.matrix, bhat=bhat,
                eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2)
            # Global test
            lrt.gl <- 2*(llkl.full - llkl_null_global(wgt.matrix=wgt.matrix, bhat=bhat, eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2))
            return.list$p.global <- 0.95*pchisq(lrt.gl,df=2,lower.tail=FALSE) +
                0.05*pchisq(lrt.gl,df=3,lower.tail=FALSE)
            # Dynamic test
            lrt.dyn <- 2*(llkl.full - llkl_null_dynamic(wgt.matrix=wgt.matrix, bhat=bhat, eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2))
            return.list$p.dynamic <- 0.95*pchisq(lrt.dyn,df=1,lower.tail=FALSE) +
                0.05*pchisq(lrt.dyn,df=2,lower.tail=FALSE)
            # Non-linear test
            lrt.nl <- 2*(llkl.full - llkl_null_nonlinear(wgt.matrix=wgt.matrix, bhat=bhat, eigR=eigR, ngwas=ngwas, omegal=omegal, omega=omega, omega2=omega2))
            return.list$p.nonlinear <- 0.05*pchisq(lrt.nl,df=1,lower.tail=FALSE)
        }

        return(return.list)
    } else{
        return(NULL)
    }

}

# Association testing using scalar-on-function regression
#' @param sumstat A data frame of GWAS summary statistics. Columns: SNP, A1, A2, Z
#' @param wgtlist A data frame of gene annotation. Five columns: ID=geneid_short, CHR=seqname, P0=start, P1=end, tss.
#' @param weights_pred A list where each each element is an output of compute_weights_flasso. The number and order of genes follow wgtlist. Named by gene ID (same as wgtlist$ID).
#' @param bim_train bim file from the training data where prediction models are built. Should include all the SNPs in wgt.
#' @param genos Reference data in plink format used to compute LD. Can be read into R using plink2R. A list of length 3: bed, bim, fam. Usually 1000 Genomes.
#' @param n_gwas GWAS sample size, can be SNP-specific or one integer.
#' @param opt Options. max_impute (maximum proportion of missing SNPs). min_r2pred
#' @param eta Tuning parameter for smoothness penalty. Default NULL, eta needs to be tuned using data thinning. Either a single number or the same length as weights_pred
sctwas_assoc_fda <- function(
        sumstat, wgtlist, weights_pred, bim_train, genos, ngwas,
        opt=list(max_impute=0.6, min_r2pred=0.8, degree_beta=3, knots_beta=seq(0.05,0.95,by=0.05),
                 logsigma2vec=seq(-20,8,by=0.2)))
{
    # Extract results from compute_weights_flasso. output.
    wgt <- lapply(weights_pred, function(x) x$Wmat[rowSums(x$Wmat!=0)>0,,drop=FALSE])

    # Load in list of weights
    chr = unique(wgtlist$CHR)

    # Sample size
    N = nrow(wgtlist) # Number of genes
    out.tbl <- wgtlist %>% mutate(sigma2=NA, p.global=NA, p.dynamic=NA, p.nonlinear=NA)

    # Match summary data to input, record NA where summary data is missing
    m = match( genos$bim[,2] , sumstat$SNP )
    sum.missing = is.na(m)
    sumstat = sumstat[m,]
    sumstat$SNP = genos$bim[,2]
    sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
    sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]

    # QC / allele-flip the input and output
    # Between the summary stats and reference genome (genos)
    qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )

    # Flip Z-scores for mismatching alleles
    sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
    sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
    sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]

    # Remove strand ambiguous SNPs (if any)
    if ( sum(!qc$keep) > 0 ) {
        genos$bim = genos$bim[qc$keep,]
        genos$bed = genos$bed[,qc$keep]
        sumstat = sumstat[qc$keep,]
    }

    FAIL.ctr = 0
    betal <- var.betal <- beta <- var.beta <- data.gene.qc <- vector("list", length=nrow(wgtlist))
    names(betal) <- names(var.betal) <- names(beta) <- names(var.beta) <- names(data.gene.qc) <- wgtlist$ID

    for (w in 1:nrow(wgtlist)){
        if (w%%50==0) print(paste("Processing",w,wgtlist$ID[w]))
        data.gene.qc.w <- qc_impute_onegene(wgt, weights_pred, wgtlist=wgtlist, w, bim_train, genos, sumstat, opt=opt)

        if (!is.null(data.gene.qc.w))
            data.gene.qc[[w]] <- data.gene.qc.w
    }

    # Association analysis for each gene
    for ( w in 1:nrow(wgtlist) ) {
        if (w%%1==0) print(paste("Fitting",w,wgtlist$ID[w]))

        if (!is.null(data.gene.qc[[w]])){
            bhat.gene <- data.gene.qc[[w]]$z/sqrt(ngwas)

            res.onegene <- assoc_onegene(
                wgt.matrix=data.gene.qc[[w]]$wgt.matrix, bhat=bhat.gene,
                R=data.gene.qc[[w]]$R, ngwas=ngwas, omegal=data.gene.qc[[w]]$omegal,
                omega=data.gene.qc[[w]]$omega, omega2=data.gene.qc[[w]]$omega2,
                method="ml", logsigma2vec=opt$logsigma2vec, test=TRUE)

            if (!is.null(res.onegene)){
                betal[[w]] <- res.onegene$betal
                var.betal[[w]] <- res.onegene$var.betal
                beta[[w]] <- res.onegene$betahat
                var.beta[[w]] <- res.onegene$var.betahat

                out.tbl$sigma2[w] <- res.onegene$sigma2
                out.tbl$p.global[w] <- res.onegene$p.global
                out.tbl$p.dynamic[w] <- res.onegene$p.dynamic
                out.tbl$p.nonlinear[w] <- res.onegene$p.nonlinear
            }

        } else{
            FAIL.ctr = FAIL.ctr + 1
        }
    }

    ind.keep <- !sapply(beta, is.null) & !sapply(var.beta, is.null) &
        !sapply(betal, is.null) & !sapply(var.betal, is.null)
    betal <- betal[ind.keep]
    var.betal <- var.betal[ind.keep]
    beta <- beta[ind.keep]
    var.beta <- var.beta[ind.keep]

    out.tbl <- out.tbl[ind.keep,]
    weights_pred.keep <- weights_pred[ind.keep]
    out.tbl$degree <- opt$degree_beta # sapply(weights_pred.keep, function(x) x$degree)
    knots <- lapply(1:nrow(out.tbl), function(x) opt$knots_beta)

    cat("Analysis completed.\n")
    cat("NOTE:",FAIL.ctr,"/",nrow(wgtlist),"genes were skipped\n")
    if ( FAIL.ctr / nrow(wgtlist) > 0.1 ) {
        cat("If a large number of genes were skipped, verify that your GWAS Z-scores, expression weights, and LDREF data use the same SNPs (or nearly)\n")
        cat("Or consider pre-imputing your summary statistics to the LDREF markers using summary-imputation software such as [https://github.com/bogdanlab/fizi]\n")
    }

    return(list(out.tbl=out.tbl, betal=betal, var.betal=var.betal, beta=beta, var.beta=var.beta, knots=knots))
}

# Test at each psuedotime point
pointwise_test <- function(sctwas_out, pt_grid=seq(0,1,by=0.01)){
    betapw <- betapw.se <-
        matrix(NA, nrow=nrow(sctwas_out$out.tbl), ncol=length(pt_grid),
               dimnames=list(sctwas_out$out.tbl$ID, paste0("pt_",1:length(pt_grid))))

    for (w in 1:nrow(sctwas_out$out.tbl)){
        ptbs <- bs(pt_grid, knots=sctwas_out$knots[[w]],
                   degree=sctwas_out$out.tbl$degree[w], intercept=TRUE, Boundary.knots=c(0,1))
        # ptbs <- cbind(1,pt_grid,ptbs[,-c(1,ncol(ptbs))])
        ptbs <- ptbs[,-c(1,ncol(ptbs))]
        betapw[w,] <- ptbs %*% sctwas_out$beta[[w]] +
            sctwas_out$betal[[w]][1] + sctwas_out$betal[[w]][2]*pt_grid

        ptlin <- cbind(1,pt_grid)
        betapw.se[w,] <- sqrt(diag(ptbs%*%sctwas_out$var.beta[[w]]%*%t(ptbs) +
                                       ptlin%*%sctwas_out$var.betal[[w]]%*%t(ptlin)))
    }

    return(list(betapw=betapw, betapw.se=betapw.se))
}

# Test for non-zero effects in an interval of pseudotime
interval_test <- function(sctwas_out, interval){
    pval <- rep(NA, nrow(sctwas_out$out.tbl))
    names(pval) <- sctwas_out$out.tbl$ID
    pt_grid <- seq(interval[1], interval[2], length.out=100)
    pt_grid <- (pt_grid[1:(length(pt_grid)-1)]+pt_grid[2:length(pt_grid)])/2
    grid_space <- mean(diff(pt_grid))

    for (w in 1:nrow(sctwas_out$out.tbl)){
        ptbs <- bs(pt_grid, knots=sctwas_out$knots[[w]],
                   degree=sctwas_out$out.tbl$degree[w], intercept=TRUE, Boundary.knots=c(0,1))
        omega_int <- t(ptbs) %*% ptbs * grid_space
        ind.nz <- which(rowSums(omega_int!=0)>0) # Indices of non-zero rows in omega_int
        omega_int_sqrt <- chol(omega_int[ind.nz,ind.nz])
        beta.nz <- sctwas_out$beta[[w]][ind.nz]
        var.beta.nz <- sctwas_out$var.beta[[w]][ind.nz,ind.nz]

        chisqstat <- t(beta.nz) %*% t(omega_int_sqrt) %*%
            solve(omega_int_sqrt%*%var.beta.nz%*%t(omega_int_sqrt),
                  omega_int_sqrt%*%beta.nz)

        pval[w] <- pchisq(chisqstat, df=length(ind.nz), lower.tail=FALSE)
    }

    return(pval)
}

## Integrated function to run scTWAS association tests
assoc_tests_fda <- function(sctwas_out, pt_grid=seq(0,1,by=0.01)){

    sctwas_out$PointwiseTest <- pointwise_test(sctwas_out=sctwas_out, pt_grid=pt_grid)
    sctwas_out$out.tbl$fdr.global <- p.adjust(sctwas_out$out.tbl$p.global, method="fdr")
    sctwas_out$out.tbl$fdr.dynamic <- p.adjust(sctwas_out$out.tbl$p.dynamic, method="fdr")
    sctwas_out$out.tbl$fdr.nonlinear <- p.adjust(sctwas_out$out.tbl$p.nonlinear, method="fdr")

    return(sctwas_out)
}
