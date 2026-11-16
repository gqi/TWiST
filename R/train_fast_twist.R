#' Train prediction model for genetically regulated gene expression using fast-TWiST
#' @description This function trains a model using single-cell expression eQTL data to predict genetically regulated gene expression (GReX) from cis-SNPs. It fits a penalized Poisson functional regression model using individual-level summaries. This is a Stage-1 model for TWiST.
#' @param y Expression count data for one gene. Integer vector of length ncells (number of cells).
#' @param geno Genotype of cis-SNPs of the gene. Numerical matrix with individuals in rows and SNPs in columns.
#' @param donor Vector of length ncells giving the row number in geno for each cell.
#' @param pt Pseudotime of the cells, scaled to rank/(number of cells) between 0 and 1.
#' @param libsize Library size of the cells. Vector of length ncells. Default is 1 for all cells.
#' @param covar Covariates to be adjusted. Numerical matrix with ncells rows, or NULL.
#' @param knots Internal knots for B-spline basis functions. Default is c(0.25, 0.5, 0.75).
#' @param degree Degree of B-spline basis functions. Default is 3 (cubic B-spline).
#' @param lambda A single penalty parameter. If NULL, selected by cross-validation.
#' @param nlambda Number of penalty values in cross-validation. Default is 20.
#' @param alpha Mixing parameter for the group elastic-net penalty. Default is 0.5; 1 gives group lasso.
#' @param nfolds Number of cross-validation folds. Default is 5.
#' @param fold Fold assignments for the cells. Vector of length ncells with values from 1 to nfolds.
#' @param eps Convergence tolerance. Default is 1e-4.
#' @param max.iter Maximum number of iterations per solver call, shared across its penalty path. Default is 10000.
#' @param ncores Number of cores for cross-validation. Default is 1.
#'
#' @return A list containing:
#' \item{Wmat}{SNP-effect coefficients. Matrix of (number of SNPs) x (number of B-spline basis functions), in the original spline basis.}
#' \item{knots, degree}{Internal knots and degree of the B-spline basis.}
#' \item{intercept}{Additional intercept fitted with the genetic effects, excluding the fixed offset.}
#' \item{lambda, alpha}{Selected penalty parameter and mixing parameter.}
#' \item{cv}{Cross-validation results, or NULL when lambda is supplied.}
#' \item{covariate_coefficients}{Coefficients of the covariate-only Poisson model, including its intercept.}
#' Other entries record the sample size, convergence, and fitting method.
#'
#' @details Covariate effects are estimated first and held fixed as an offset, including log library size. SNP effects are represented using B-splines and fitted with a group elastic-net penalty. Individual-level spline summaries are used without constructing the full cell-by-SNP spline design matrix.
#' @importFrom stats glm.fit poisson glm residuals coef pcauchy
#' @importFrom parallel mclapply
#' @useDynLib TWiST, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
fast_twist_train_model <- function(y, geno, donor, pt, libsize=rep(1,length(y)),
        covar=NULL, knots=c(0.25,0.5,0.75), degree=3,
        lambda=NULL, nlambda=20, alpha=0.5, nfolds=5,
        fold=sample(rep(1:nfolds, length.out=length(y))),
        eps=1e-4, max.iter=10000, ncores=1){
    n <- length(y)

    # Fit covariates first and use the linear predictor as offset
    X <- cbind("(Intercept)"=rep(1,n), covar)
    covariate_fit <- glm.fit(X, y, family=poisson(), offset=log(libsize))
    offset <- covariate_fit$linear.predictors

    # B-spline bases and standardization of SNP groups
    B <- bs(pt, knots=knots, degree=degree, intercept=TRUE, Boundary.knots=c(0,1))
    K <- ncol(geno)
    M <- ncol(B)
    full <- fast_twist_standardize(geno, B, donor)

    # Cross validation to select optimal lambda parameter
    cv <- NULL
    if (is.null(lambda)){
        # lambda sequence starts from the max gradient under the null model, same as grpreg
        null_fit <- glm(y ~ 1, family="poisson", offset=offset)
        r <- residuals(null_fit, "working") * null_fit$weights
        H <- matrix(0, nrow(geno), M)
        H[sort(unique(donor)), ] <- rowsum(B * r, donor)
        GH <- crossprod(geno, H)
        rsum <- sum(r)
        z <- sapply(1:K, function(k){
            ix <- (k - 1) * M + 1:M
            grad <- crossprod(full$T_array[, , k], (GH[k, ] - full$center[ix] * rsum) / full$scale[ix])
            sqrt(sum(grad^2))
        })
        zmax <- max((z / full$m)[full$m > 0])
        lambda.max <- zmax / n / alpha
        lambda.min <- if (n > K * M) 1e-4 else 0.05
        lambda <- exp(seq(log(lambda.max), log(lambda.min * lambda.max), length=nlambda))

        fits <- mclapply(1:nfolds, function(i){
            train <- fold != i
            test <- !train
            representation <- fast_twist_standardize(geno, B[train, , drop=FALSE],
                donor[train], full$feature_off, full$Q_array)
            fit <- flasso_compressed_rep_cd_fit_cpp(
                G=geno, B=B[train, , drop=FALSE], donor_index=donor[train],
                y=y[train], lambda=lambda, feature_off=representation$feature_off,
                Q_array=representation$Q_array, group_multiplier=full$m,
                alpha=alpha, eps=eps, max_iter=max.iter, offset_vec=offset[train])

            # Poisson deviance of the held-out cells for each lambda
            mu <- matrix(NA, sum(test), length(lambda))
            for (l in 1:length(lambda)){
                beta <- fit$beta_ortho[, l]
                intercept <- fit$intercept_ortho[l] + sum(representation$feature_off * beta)
                W <- matrix(0, K, M)
                for (k in 1:K){
                    ix <- (k - 1) * M + 1:M
                    W[k, ] <- drop(representation$Q_array[, , k] %*% beta[ix])
                }
                F <- geno %*% W
                eta <- offset[test] + intercept +
                    rowSums(B[test, , drop=FALSE] * F[donor[test], , drop=FALSE])
                mu[, l] <- exp(eta)
            }
            err <- 2 * (mu - y[test])
            nz <- y[test] > 0
            err[nz, ] <- err[nz, , drop=FALSE] +
                2 * y[test][nz] * log(y[test][nz] / mu[nz, , drop=FALSE])
            err[, !fit$converged] <- NA
            list(sum=colSums(err), sumsq=colSums(err^2),
                 converged=fit$converged, iter=fit$iter)
        }, mc.cores=min(ncores, nfolds), mc.set.seed=FALSE)
        cve <- Reduce(`+`, lapply(fits, `[[`, "sum")) / n
        sumsq <- Reduce(`+`, lapply(fits, `[[`, "sumsq"))
        cvse <- sqrt(pmax((sumsq - n * cve^2) / (n - 1), 0) / n)
        valid <- is.finite(cve)
        if (!any(valid)) stop("No regularization value has valid converged fits in all CV folds.")
        selected <- which(valid)[which.min(cve[valid])]
        cv <- list(lambda=lambda, cve=cve, cvse=cvse, fold=fold,
                   selected=selected, converged=lapply(fits, `[[`, "converged"),
                   iter=lapply(fits, `[[`, "iter"))
        lambda <- cv$lambda[cv$selected]
    }

    # Refit using all cells
    fit <- flasso_compressed_rep_cd_fit_cpp(
        G=geno, B=B, donor_index=donor, y=y, lambda=lambda,
        feature_off=full$feature_off, Q_array=full$Q_array,
        group_multiplier=full$m, alpha=alpha, eps=eps,
        max_iter=max.iter, offset_vec=offset)
    if (!fit$converged[1] || any(!is.finite(fit$beta_ortho)))
        stop("fast-TWiST fit failed: nonconvergence or non-finite coefficients.")
    # Transform coefficients back to the original B-spline bases
    Wmat <- matrix(0, K, M, dimnames=list(colnames(geno), paste0("bs",1:M)))
    for (k in 1:K){
        ix <- (k - 1) * M + 1:M
        Wmat[k, ] <- drop(full$Q_array[, , k] %*% fit$beta_ortho[ix, 1])
    }
    intercept <- fit$intercept_ortho[1] + sum(full$feature_off * fit$beta_ortho[, 1])
    list(Wmat=Wmat, knots=knots, degree=degree, intercept=intercept,
         lambda=lambda, alpha=alpha, cv=cv, iter=fit$iter[1],
         converged=fit$converged[1], n=length(y),
         covariate_coefficients=covariate_fit$coefficients, method="fast-TWiST")
}

#' Standardize SNP groups using individual-level spline summaries
#' @param geno_donor Genotype matrix with individuals in rows and SNPs in columns.
#' @param ptbs B-spline basis matrix with cells in rows.
#' @param donor_index Genotype row number for each cell.
#' @param feature_off Feature offsets from an earlier standardization, or NULL.
#' @param Q_array Basis transformations from an earlier standardization, or NULL.
#' @param svd_tol Tolerance for detecting rank-deficient SNP groups.
#' @return A list of feature offsets, basis transformations, centering and scaling values, and group penalty multipliers.
fast_twist_standardize <- function(geno_donor, ptbs, donor_index,
        feature_off=NULL, Q_array=NULL, svd_tol=1e-10){
    n_snp <- ncol(geno_donor)
    n_basis <- ncol(ptbs)
    n <- nrow(ptbs)
    S <- matrix(0, nrow(geno_donor), n_basis)
    M <- array(0, dim=c(nrow(geno_donor), n_basis, n_basis))
    ids <- sort(unique(donor_index))
    S[ids, ] <- rowsum(ptbs, donor_index)
    for (k in 1:n_basis)
        M[ids, k, ] <- rowsum(ptbs * ptbs[, k], donor_index)
    if (is.null(feature_off)){
        feature_off <- numeric(n_snp * n_basis)
        Q_array <- array(0, dim=c(n_basis, n_basis, n_snp))
        for (k in 1:n_snp) Q_array[, , k] <- diag(n_basis)
    }

    group_multiplier <- numeric(n_snp)
    out_off <- numeric(n_snp * n_basis)
    out_Q <- array(0, dim=c(n_basis, n_basis, n_snp))
    center <- numeric(n_snp * n_basis)
    scale <- rep(1, n_snp * n_basis)
    T_array <- array(0, dim=c(n_basis, n_basis, n_snp))

    for (j in 1:n_snp){
        idx <- ((j - 1) * n_basis + 1):(j * n_basis)
        gj <- geno_donor[, j]
        off0 <- feature_off[idx]
        Q0 <- Q_array[, , j]
        z_sum_basis <- drop(crossprod(gj, S))
        z_cross_basis <- matrix(0, n_basis, n_basis)
        gj2 <- gj * gj
        for (k in 1:n_basis){
            for (l in 1:n_basis){
                z_cross_basis[k, l] <- sum(gj2 * M[, k, l])
            }
        }
        z_sum <- drop(crossprod(z_sum_basis, Q0))
        z_cross <- crossprod(Q0, z_cross_basis %*% Q0)
        raw_cross <- n * tcrossprod(off0) + tcrossprod(off0, z_sum) +
            tcrossprod(z_sum, off0) + z_cross
        mu <- off0 + z_sum / n
        sc <- sqrt(pmax(diag(raw_cross) / n - mu^2, 0))
        gram <- raw_cross - n * tcrossprod(mu)
        denom <- pmax(sc, 1e-6)
        gram <- sweep(sweep(gram, 1, denom, "/"), 2, denom, "/") / n
        eig <- eigen((gram + t(gram)) / 2, symmetric=TRUE)
        if (any(sc <= 1e-6) ||
            min(eig$values) <= sqrt(.Machine$double.eps) * max(eig$values)){
            # Rank deficient group: use SVD of the cell-level block, same as grpreg
            Z <- sweep((ptbs * gj[donor_index]) %*% Q0, 2, off0, "+")
            mu <- colMeans(Z)
            Z <- sweep(Z, 2, mu, "-")
            sc <- sqrt(colSums(Z^2) / n)
            keep <- which(sc > 1e-6)
            if (!length(keep)) next
            SVD <- svd(sweep(Z[, keep, drop=FALSE], 2, sc[keep], "/"), nu=0)
            r <- which(SVD$d > svd_tol)
            Tj <- matrix(0, n_basis, n_basis)
            Tj[keep, seq_along(r)] <- sweep(SVD$v[, r, drop=FALSE],
                                           2, sqrt(n) / SVD$d[r], "*")
            sc[-keep] <- 1
            group_multiplier[j] <- sqrt(length(r))
        } else{
            V <- eig$vectors
            for (k in 1:ncol(V)){
                imax <- which.max(abs(V[, k]))
                if (V[imax, k] < 0) V[, k] <- -V[, k]
            }
            Tj <- sweep(V, 2, 1 / sqrt(eig$values), "*")
            group_multiplier[j] <- sqrt(n_basis)
        }
        A <- sweep(Tj, 1, sc, "/")
        out_Q[, , j] <- Q0 %*% A
        out_off[idx] <- drop((off0 - mu) %*% A)
        center[idx] <- mu
        scale[idx] <- sc
        T_array[, , j] <- Tj
    }
    list(feature_off=out_off, Q_array=out_Q,
         center=center, scale=scale, T_array=T_array,
         m=group_multiplier)
}
