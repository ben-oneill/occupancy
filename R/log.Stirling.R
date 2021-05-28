#' Logarithms of the Stirling numbers of the second kind
#'
#' \code{log.Stirling} returns a matrix of the logarithms of the Stirling numbers of the second kind.
#'
#' This function computes a matrix of the logarithms of the Stirling numbers of the second kind.  The function
#' allows the user to give a non-centrality parameter for the non-central Stirling numbers.
#'
#' @usage \code{log.Stirling(n, k, ncp = 0)}
#' @param n A vector of non-negative integer values
#' @param k A vector of non-negative integer values
#' @param ncp Non-centrality parameter (non-negative numeric value)
#' @return If all inputs are correctly specified then the output will be a matrix containing the logarithms of
#' the Stirling numbers of the second kind

log.Stirling <- function(n, k, ncp = 0) {

  #Check that argument values are single numeric values
  if (!is.numeric(n))                        stop('Error: Argument n is not numeric')
  if (!is.numeric(k))                        stop('Error: Argument k is not numeric')
  if (!is.numeric(ncp))                      stop('Error: Argument ncp is not numeric')
  if (length(ncp) != 1)                      stop('Error: Argument ncp should be a single nonnegative number')

  #Set parameters
  N <- as.integer(n)
  K <- as.integer(k)

  #Check that arguments are in allowable range
  if (any(n != N))                           stop('Error: Values in n must be non-negative integers')
  if (any(k != K))                           stop('Error: Values in k must be non-negative integers')
  if (min(N) < 0)                            stop('Error: Values in n must be non-negative integers')
  if (min(K) < 0)                            stop('Error: Values in k must be non-negative integers')
  if (ncp < 0)                               stop('Error: Argument ncp must be non-negative')

  #Set parameters
  n <- max(N)
  k <- max(K)

  #Set log-Stirling matrix
  LOGSTIRLING <- matrix(-Inf, nrow = n+1, ncol = k+1)
  rownames(LOGSTIRLING) <- sprintf('n[%s]', 0:n)
  colnames(LOGSTIRLING) <- sprintf('k[%s]', 0:k)

  #Compute Stirling numbers of the second kind
  #Generate base log-Stirling numbers
  LOGSTIRLING[1,1] <- 0
  if ((ncp > 0)&(n > 0)) {
  for (nn in 1:n) {
    LOGSTIRLING[nn+1, 1] <- nn*log(ncp) } }

  #Generate log-Stirling numbers via recursion
  LOGSTIRLING[1,1] <- 0
  for (nn in 1:n) {
  for (kk in 1:min(k,nn)) {
    T1 <- log(kk + ncp) + LOGSTIRLING[nn, kk+1]
    T2 <- LOGSTIRLING[nn, kk]
    LOGSTIRLING[nn+1, kk+1] <- matrixStats::logSumExp(c(T1, T2)) } }

  #Set output
  OUT <- LOGSTIRLING[N+1, K+1, drop = FALSE]
  attr(OUT, 'Description') <- 'Log-Stirling numbers of the second kind'
  attr(OUT, 'ncp') <- ncp

  #Return output
  OUT }
