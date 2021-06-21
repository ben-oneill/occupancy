#' occupancy: A package for computing with occupancy distributions.
#'
#' The occupancy package provides standard mass functions, distribution functions, quantile functions
#' and random generation for various distributions.  These distributions occur in cases where balls are
#' randomly allocated to bins.
#'
#'
#' @docType package
#' @name occupancy
NULL




#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	logical; if TRUE (default), probabilities are \eqn{P[X ≤ x]} otherwise, \eqn{P[X > x]}.
#' @keywords internal
.inheritparams <- function(x, q, p, n, log, log.p, lower.tail) NULL
