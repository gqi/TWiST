#' Train prediction model for genetically regulated gene expression using FPCA
#' @description This function predicts genetically regulated gene expression (GReX) from cis-SNPs using functional principal component analysis (FPCA). Individual expression trajectories are estimated using Poisson mixed-effects models within pseudotime bins. Component scores are predicted from genotypes using elastic net regression. This is a Stage-1 model for TWiST.
#' @param y Expression count data for one gene. Integer vector of length ncells (number of cells).
#' @param geno Genotype of cis-SNPs of the gene. Numerical matrix with individuals in rows and SNPs in columns.
#' @param donor Vector of length ncells giving the row number in geno for each cell.
#' @param pt Pseudotime of the cells, scaled to rank/(number of cells) between 0 and 1.
#' @param libsize Library size of the cells. Vector of length ncells. Default is 1 for all cells.
#' @param covar Covariates to be adjusted. Numerical matrix with ncells rows, or NULL.
#' @param nbins Number of approximately equal-cell-count pseudotime bins before merging zero-count bins. Default is 50.
#' @param pve Proportion of trajectory variation explained by the retained components. Default is 0.99.
#' @param alpha Elastic-net mixing parameter for score regression. Default is 0.1.
#' @param nfolds Number of cross-validation folds for score regression. Default is 10.
#' @param fold Fold assignments for individuals, ordered by their observed genotype row numbers. The same folds are used for all components.
#'
#' @return A list containing:
#' \item{Wmat}{Genotype coefficients for the component scores. Matrix of (number of SNPs) x (number of retained components).}
#' \item{efunctions}{Eigenfunctions evaluated at the bin coordinates, with bins in rows and components in columns. Together with Wmat, these define the SNP-effect curves.}
#' \item{argvals, bin_breaks}{Pseudotime coordinates and boundaries of the bins.}
#' \item{mean, intercepts_pc}{Mean expression trajectory and score-regression intercepts.}
#' \item{lambda, alpha}{Selected penalty parameters and elastic-net mixing parameter.}
#' \item{fold, donor_rows}{Cross-validation folds and corresponding genotype row numbers.}
#' \item{covariate_coefficients}{Coefficients of the covariate-only Poisson model, including its intercept.}
#' Other entries record the number of components, binning and FPCA settings, sample size, and fitting method.
#'
#' @details Covariate effects are estimated first and held fixed as an offset, including log library size. A Poisson mixed-effects model with an individual-specific random intercept is fitted within each bin. The random effects across bins define the individual expression trajectories. FPCA is fitted using fpca.sc with its default smoothing settings. Component scores are standardized before elastic net regression, and the fitted coefficients are returned on the original score scale.
#' @importFrom lme4 glmer ranef
#' @importFrom refund fpca.sc
#' @importFrom glmnet cv.glmnet
#' @importFrom stats median sd
#' @export
fpca_train_model <- function(y, geno, donor, pt, libsize=rep(1,length(y)),
        covar=NULL, nbins=50, pve=0.99, alpha=0.1, nfolds=10,
        fold=sample(rep(1:nfolds, length.out=length(unique(donor))))){
    n <- length(y)

    # Fit covariates first and use the linear predictor as offset
    X <- cbind("(Intercept)"=rep(1,n), covar)
    covariate_fit <- glm.fit(X, y, family=poisson(), offset=log(libsize))
    offset <- covariate_fit$linear.predictors

    # Pseudotime bins with equal numbers of cells
    bin <- ntile(pt, nbins)
    # Bins with zero counts are merged into the preceding nonzero bin (leading ones into the first)
    nonzero <- which(tapply(y, bin, sum) > 0)
    bin <- pmax(1, findInterval(bin, nonzero))
    nbin <- max(bin)
    ids <- sort(unique(donor))
    curves <- matrix(NA, length(ids), nbin, dimnames=list(ids, NULL))
    argvals <- sapply(1:nbin, function(j) median(pt[bin == j]))
    ends <- sapply(1:nbin, function(j) max(pt[bin == j]))
    breaks <- c(0, ends[-nbin], 1)
    # Poisson mixed model in each bin
    for (j in 1:nbin){
        ix <- bin == j
        df <- data.frame(y=y[ix], id=factor(donor[ix]), offset=offset[ix])
        fit <- glmer(y ~ 1 + (1 | id), data=df, family=poisson(), offset=offset, nAGQ=0)
        effects <- ranef(fit)$id
        curves[match(rownames(effects), rownames(curves)), j] <- effects[, 1]
    }

    # FPCA, then elastic net for the score of each component
    fp <- fpca.sc(Y=curves, argvals=argvals, pve=pve)
    scores <- fp$scores
    G <- geno[ids, , drop=FALSE]
    L <- ncol(scores)
    theta <- matrix(0, ncol(G), L, dimnames=list(colnames(G), paste0("pc",1:L)))
    intercepts <- lambda <- numeric(L)
    for (l in 1:L){
        center <- mean(scores[, l])
        scale <- sd(scores[, l])
        z <- (scores[, l] - center) / scale
        cv <- cv.glmnet(G, z, family="gaussian", alpha=alpha,
                        nfolds=nfolds, foldid=fold)
        lambda[l] <- cv$lambda.min
        coef <- as.numeric(coef(cv, s="lambda.min"))
        theta[, l] <- coef[-1] * scale
        intercepts[l] <- coef[1] * scale + center
    }

    list(Wmat=theta, efunctions=fp$efunctions, argvals=argvals,
         bin_breaks=breaks, mean=fp$mu, intercepts_pc=intercepts,
         lambda=lambda, alpha=alpha, fold=fold, donor_rows=ids,
         npc=L, pve=pve, nbins=nbin, n=length(y),
         covariate_coefficients=covariate_fit$coefficients, method="FPCA")
}
