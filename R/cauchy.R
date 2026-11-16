#' Combine fast-TWiST and FPCA P values
#' @description This function combines paired P values using the Cauchy combination test with equal weights. Global and dynamic tests are combined separately.
#' @param p_fast Vector of fast-TWiST P values.
#' @param p_fpca Vector of FPCA P values for the same genes and hypothesis, in the same order.
#' @return A vector of combined P values. Entries with a missing or non-finite component P value are returned as NA.
#' @export
cauchy_combine <- function(p_fast, p_fpca){
    p <- cbind(p_fast, p_fpca)
    out <- rep(NA, nrow(p))
    valid <- rowSums(is.finite(p)) == 2
    p <- p[valid, , drop=FALSE]
    p <- pmin(pmax(p, 1e-300), 1 - 1e-16)
    z <- tan((0.5 - p) * pi)
    small <- p < 1e-15
    z[small] <- 1 / (pi * p[small])
    stat <- rowMeans(z)
    ans <- pcauchy(stat, lower.tail=FALSE)
    large <- stat > 1e15
    ans[large] <- 1 / (pi * stat[large])
    out[valid] <- ans
    out
}
