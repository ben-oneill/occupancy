#' Cumulative distribution function of the extended occupancy distribution
#'
#' \code{pocc} returns the probability or log-probability values for the arguments.
#'
#' This function computes probabilities or log-probabilities from the cumulative distribution function of the
#' extended occupancy distribution.  Further details on the distribution can be found in the following paper:
#'
#' O'Neill, B. (2021) Three distributions in the extended occupancy problem.
#'
#' @usage \code{pocc(x, size, space, prob, approx = FALSE, log.p = FALSE, lower.tail = TRUE)}
#' @param x A vector of numeric values to be used as arguments for the mass function
#' @param size The size parameter for the occupancy distribution (number of balls)
#' @param space The space pararmeter for the occupancy distribution (number of bins)
#' @param prob The probability parameter for the occupancy distribution (probability of ball occupying its bin)
#' @param approx A logical value specifying whether to use an approximation to the distribution
#' @param log.p A logical value specifying whether results should be returned as log-probabilities
#' @param lower.tail A logical value specifying whether results are from the cumulative distribution function
#' or the corresponding survival function
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are integers)
#' then the output will be a vector of probabilities/log-probabilities corresponding to the vector argument x

pocc <- function(x, size, space, prob = 1, approx = FALSE, log.p = FALSE, lower.tail = TRUE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(x))                                   stop('Error: Argument x is not numeric')
  if (!is.numeric(size))                                stop('Error: Size parameter is not numeric')
  if (!is.numeric(space))                               stop('Error: Space parameter is not numeric')
  if (!is.numeric(prob))                                stop('Error: Probability parameter is not numeric')
  if (!is.logical(approx))                              stop('Error: approx option is not a logical value')
  if (!is.logical(log.p))                               stop('Error: log.p option is not a logical value')
  if (!is.logical(lower.tail))                          stop('Error: lower.tail option is not a logical value')

  #Check that parameters are atomic
  if (length(size)  != 1)                               stop('Error: Size parameter should be a single number')
  if (length(space) != 1)                               stop('Error: Space parameter should be a single number')
  if (length(prob)  != 1)                               stop('Error: Probability parameter should be a single number')
  if (length(approx)   != 1)                            stop('Error: approx option should be a single logical value')
  if (length(log.p) != 1)                               stop('Error: log.p option should be a single logical value')
  if (length(lower.tail) != 1)                          stop('Error: lower.tail option should be a single logical value')

  #Set parameters
  n   <- as.integer(size)
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }
  MAX <- min(n,m)

  #Check that parameters are in allowable range
  if (size != n)                                        stop('Error: Size parameter is not an integer')
  if (n < 0)                                            stop('Error: Size parameter must be non-negative')
  if (space != m)                                       stop('Error: Space parameter is not an integer')
  if (m <= 0)                                           stop('Error: Space parameter must be positive')
  if ((prob < 0)|(prob > 1))                            stop('Error: Probability parameter is not between zero and one')

  #Create output vector
  MAX <- min(n, m)
  CUMOCC <- rep(-Inf, length(x))

  #Compute for trivial case where n = 0 or prob = 0
  if ((n == 0)|(prob == 0)) {
    IND <- (x >= 0)
    CUMOCC[IND] <- 0
    if (log) { return(CUMOCC) } else { return(exp(CUMOCC)) } }

  #Compute for trivial case where m = Inf and prob > 0
  if (m == Inf) {
    CUMOCC <- pbinom(x, size = n, prob = prob, log.p = TRUE)
    if (log) { return(CUMOCC) } else { return(exp(CUMOCC)) } }

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
    LOGS <- rep(-Inf, MAX+1)
    for (k in 0:MAX) {
      LOGS[k+1] <- n*log(prob) - n*log(m) + lchoose(m,k) + lfactorial(k) + LOGSTIRLING[n+1, k+1] }
    LOGS <- LOGS - matrixStats::logSumExp(LOGS) }

  if (approx) {

    #Compute normal approximation to the occupancy distribution
    E1   <- (1 - prob/m)^n
    E2   <- (1 - 2*prob/m)^n
    MEAN <- m*(1 - E1)
    VAR  <- m*((m-1)*E2 + E1 - m*E1^2)
    LOGS <- dnorm(0:MAX, mean = MEAN, sd = sqrt(VAR), log = TRUE)
    LOGS <- LOGS - matrixStats::logSumExp(LOGS) }

  #Generate the log-probabilities for the cumulative distribution
  CUMLOGS <- rep(-Inf, MAX+1)
  CUMLOGS[1] <- LOGS[1]
  for (k in 1:MAX) {
    CUMLOGS[k+1] <- matrixStats::logSumExp(c(CUMLOGS[k], LOGS[k+1])) }

  #Generate output vector
  for (i in 1:length(x)) {
    xx <- floor(x[i])
    if ((xx >= 0)&(xx <= MAX)) {
      CUMOCC[i] <- CUMLOGS[xx+1] }
    if (xx > MAX) {
      CUMOCC[i] <- 0 } }
  if (!lower.tail) { CUMOCC <- VGAM::log1mexp(-CUMOCC) }

  #Return output
  if (log.p) { CUMOCC } else { exp(CUMOCC) } }
