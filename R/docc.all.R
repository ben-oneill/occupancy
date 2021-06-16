#' Mass function of the extended occupancy distribution
#'
#' \code{docc.all} returns the probability or log-probability values for the arguments.
#'
#' This function computes probabilities or log-probabilities from the mass function of the extended occupancy
#' distribution.  Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2021) Three distributions in the extended occupancy problem.
#'
#' @param max.size The maximum size parameter for the occupancy distribution (number of balls)
#' @param space The space pararmeter for the occupancy distribution (number of bins)
#' @param prob The probability parameter for the occupancy distribution (probability of ball occupying its bin)
#' @param approx A logical value specifying whether to use the normal approximation to the occupancy distribution
#' @param log A logical value specifying whether results should be returned as log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are integers)
#' then the output will be a vector of probabilities/log-probabilities corresponding to the vector argument x
#' @rdname docc
docc.all <- function(max.size, space, prob = 1, approx = FALSE, log = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(max.size))                stop('Error: Maximum size parameter is not numeric')
  if (!is.numeric(space))                   stop('Error: Space parameter is not numeric')
  if (!is.numeric(prob))                    stop('Error: Probability parameter is not numeric')
  if (!is.logical(approx))                  stop('Error: approx option is not a logical value')
  if (!is.logical(log))                     stop('Error: log option is not a logical value')

  #Check that parameters are atomic
  if (length(max.size)  != 1)               stop('Error: Maximum size parameter should be a single number')
  if (length(space) != 1)                   stop('Error: Space parameter should be a single number')
  if (length(prob)  != 1)                   stop('Error: Probability parameter should be a single number')
  if (length(approx)   != 1)                stop('Error: approx option should be a single logical value')
  if (length(log)   != 1)                   stop('Error: log option should be a single logical value')

  #Set parameters
  n   <- as.integer(max.size)
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }
  MAX <- min(n,m)

  #Check that parameters are in allowable range
  if (max.size != n)                        stop('Error: Maximum size parameter is not an integer')
  if (n < 0)                                stop('Error: Maximum size parameter must be non-negative')
  if (space != m)                           stop('Error: Space parameter is not an integer')
  if (m <= 0)                               stop('Error: Space parameter must be positive')
  if ((prob < 0)|(prob > 1))                stop('Error: Probability parameter is not between zero and one')

  #Deal with trivial case where n = 0
  if (n == 0) {
    OCC <- matrix(-Inf, nrow = MAX+1, ncol = 1)
    rownames(OCC) <- sprintf('x[%s]', 0:MAX)
    colnames(OCC) <- 'n[0]'
    OCC[1,1] <- 0
    if (log) { return(OCC) } else { return(exp(OCC)) } }

  #Create output vector
  MAX <- min(n, m)
  OCC <- matrix(-Inf, nrow = MAX+1, ncol = n+1)
  rownames(OCC) <- sprintf('x[%s]', 0:MAX)
  colnames(OCC) <- sprintf('n[%s]', 0:n)

  #Compute for trivial case where prob = 0
  if (prob == 0) {
    OCC[1, ] <- 0
    if (log) { return(OCC) } else { return(exp(OCC)) } }

  #Compute for trivial case where m = Inf and prob > 0
  if (m == Inf) {
    for (nn in 1:n) {
      OCC[, nn] <- dbinom(0:n, size = n, prob = prob, log = TRUE) }
    if (log) { return(OCC) } else { return(exp(OCC)) } }

  #Compute for non-trivial case where m < Inf and prob > 0
  if (!approx) {

    #Compute log-probablities using recursion
    SCALE <- m*(1-prob)/prob

    #Set log-Stirling matrix and generate first row
    LOGSTIRLING <- matrix(-Inf, nrow = n+1, ncol = MAX+1)
    LOGSTIRLING[1,1] <- 0
    if ((SCALE > 0)&(n > 0)) {
      for (nn in 1:n) {
        LOGSTIRLING[nn+1, 1] <- nn*log(SCALE) } }

    #Generate subsequent rows
    for (nn in 1:n) {
      for (kk in 1:MAX) {
        T1 <- log(kk + SCALE) + LOGSTIRLING[nn, kk+1]
        T2 <- LOGSTIRLING[nn, kk]
        LOGSTIRLING[nn+1, kk+1] <- matrixStats::logSumExp(c(T1, T2)) } }

    #Generate the log-probabilities for the occupancy distribution
    for (nn in 1:n) {
    for (kk in 0:min(nn, m)) {
      OCC[kk+1, nn+1] <- nn*log(prob) - nn*log(m) + lchoose(m,kk) + lfactorial(kk) + LOGSTIRLING[nn+1, kk+1] }
      OCC[, nn+1] <- OCC[, nn+1] - matrixStats::logSumExp(OCC[, nn+1]) } }

  if (approx) {

    #Compute normal approximation to the occupancy distribution
    for (nn in 1:n) {
      E1   <- (1 - prob/m)^nn
      E2   <- (1 - 2*prob/m)^nn
      MEAN <- m*(1 - E1)
      VAR  <- m*((m-1)*E2 + E1 - m*E1^2)
      OCC[1:(min(nn, m)+1), nn+1] <- dnorm(0:min(nn, m), mean = MEAN, sd = sqrt(VAR), log = TRUE)
      OCC[, nn+1] <- OCC[, nn+1] - matrixStats::logSumExp(OCC[, nn+1]) } }

  #Return output
  if (log) { OCC } else { exp(OCC) } }
