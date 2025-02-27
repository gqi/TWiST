library(splines)
library(grpreg)
#' @param y Expression count for one gene. Vector of length ncells (number of cells).
#' @param geno_cell Genotype of cis-SNPs of the gene. Matrix with dimensions ncells x (number of SNPs). Each row is a cell and each column is a SNP.
#' @param pt Pseudotime of the cells. Pseudotime values learned from packages such as Slingshot, TSCAN should be scaled to rank/(number of cells) such that values are approximately uniformly distributed between 0 and 1.
#' @param knots Internal knots for B-spline basis functions that model SNP effects on gene expression. Not including 0 and 1.
#' @param degree Degree of B-spline basis functions, default to 3 (cubic B-spline).
#' @param lambda Tuning parameter for group lasso penalty. If null, `lambda` is selected by cross-validation.
#' @param nlambda Number of `lambda` values in cross-validation - default is 50.
#' @param libsize Library size of the cells. Vector of length ncells.
#' @param covar Covariate matrix to be adjusted. Not penalized.
compute_weights_poisson <- function(y, geno_cell, pt, knots=c(0.25,0.5,0.75), degree, lambda=NULL, nlambda=50, libsize=NULL, covar=NULL){

    # intercept needs to be set to TRUE to get the complete set of bases (the bases sum to 1 at each pseudotime point)
    # If intercept=FALSE, bs() will remove one basis
    ptbs <- bs(pt, knots=knots, degree=degree, intercept=TRUE, Boundary.knots=c(0,1))
    desmat <- generate_desmat(geno.temp=geno_cell, ptbs=ptbs)
    grp <- rep(1:ncol(geno_cell),each=ncol(ptbs))

    # Cross validation to select optimal lambda parameter
    if (is.null(lambda)){
        gr_cv <- cv.grpreg(X=cbind(log(libsize),covar,desmat), y=y,
                           group=c(rep(0,1+ncol(covar)), grp), penalty="grLasso",
                           family="poisson", nlambda=nlambda, alpha=0.5, nfolds=5)
        lambda.min <- gr_cv$lambda.min
    } else{
        lambda.min <- lambda
    }

    res.opt <- grpreg(X=cbind(log(libsize),covar,desmat), y=y,
                      group=c(rep(0,1+ncol(covar)), grp), penalty="grLasso",
                      family="poisson", lambda=lambda.min, alpha=0.5)
    # Different from glmnet, res.opt$beta from grpnet includes intercept
    res.opt$betamat <- matrix(res.opt$beta[-(1:(ncol(covar)+2))], nrow=ncol(geno_cell), byrow=TRUE,
                              dimnames=list(colnames(geno_cell),paste0("bs",1:ncol(ptbs))))
    res.opt$knots <- knots
    res.opt$degree <- degree

    if (is.null(lambda)){
        res.opt$lambda.seq <- gr_cv$lambda
        res.opt$cvm <- gr_cv$cvm
    }

    return(res.opt)
}

# Generate design matrix: geno.temp should have the same number of rows as ptbs
generate_desmat <- function(geno.temp, ptbs){
    if (nrow(geno.temp)!=nrow(ptbs)) stop("geno.temp and ptbs should have the same number of rows.")

    desmat <- matrix(NA, nrow=nrow(geno.temp), ncol=ncol(ptbs)*ncol(geno.temp))
    for (i in 1:ncol(geno.temp)){
        desmat[,(ncol(ptbs)*(i-1)+1):(ncol(ptbs)*i)] <- ptbs*geno.temp[,i]
    }

    return(desmat)
}
