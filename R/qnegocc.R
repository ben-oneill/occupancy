#' Quantile function of the negative occupancy distribution
#'
#' \code{qnegocc} returns the quantile function for the distribution.
#'
#' This function computes values from the quantile function of the negative occupancy distribution.  The computation
#' method uses a recursive algorithm to compute the cumulative probabilities and then converts this to quantiles using
#' the probabilities/log-probabilities in the input argument.  The recursive method is from the following paper:
#' in the following paper:
#'
#' O'Neill, B. (2021) An examination of the negative-occupancy distribution and the coupon-collector distribution.
#'
#' @usage \code{qnegocc(p, space, occupancy, prob, approx = FALSE, log.p = FALSE, lower.tail = TRUE)}
#' @param p A vector of numeric probability/log-probability values to be used as arguments for the quantile function
#' @param space The space parameter for the negative occupancy distribution (number of bins)
#' @param occupancy The occupancy parameter for the negative occupancy distribution (number of occupied bins)
#' @param prob The probability parameter for the negative occupancy distribution (probability of ball occupying its bin)
#' @param approx A logical value specifying whether to use an approximation for the distribution
#' @param log.p A logical value specifying whether input arguments are log-probabilities
#' @param lower.tail A logical value specifying whether probabilities are from the cumulative distribution function
#' or the corresponding survival function
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are integers)
#' then the output will be a vector of probabilities/log-probabilities corresponding to the vector argument p

qnegocc <- function(p, space = 1, occupancy = space, prob = 1, approx = FALSE, log.p = FALSE, lower.tail = TRUE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(p))                        stop('Error: Argument p is not numeric')
  if (!is.numeric(space))                    stop('Error: Space parameter is not numeric')
  if (!is.numeric(occupancy))                stop('Error: Occupancy parameter is not numeric')
  if (!is.numeric(prob))                     stop('Error: Probability parameter is not numeric')
  if (!is.logical(log.p))                    stop('Error: log.p option is not a logical value')
  if (!is.logical(lower.tail))               stop('Error: lower.tail option is not a logical value')

  #Check that parameters are atomic
  if (length(space)  != 1)                   stop('Error: Space parameter should be a single number')
  if (length(occupancy) != 1)                stop('Error: Occupancy parameter should be a single number')
  if (length(prob)  != 1)                    stop('Error: Probability parameter should be a single number')
  if (length(log.p) != 1)                    stop('Error: log.p option should be a single logical value')
  if (length(lower.tail) != 1)               stop('Error: lower.tail option should be a single logical value')

  #Set parameters
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }
  k <- as.integer(occupancy)

  #Check that parameters are in allowable range
  if (space != m)                            stop('Error: Space parameter is not an integer')
  if (m <= 0)                                stop('Error: Space parameter must be positive')
  if (occupancy != k)                        stop('Error: Occupancy parameter is not an integer')
  if (k > m)                                 stop('Error: Occupancy parameter is larger than space parameter')
  if (prob < 0)                              stop('Error: Probability parameter must be between zero and one')
  if (prob > 1)                              stop('Error: Probability parameter must be between zero and one')

  #Check that argument values are in allowable range
  if (!log.p) {
    if (min(p) < 0)             stop('Error: Probability values in p must be between zero and one')
    if (max(p) > 1)             stop('Error: Probability values in p must be between zero and one') }
  if (log.p) {
    if (max(p) > 0)             stop('Error: Log-probability values in p must be less than or equal to zero') }

  #Set maximum log-probability for quantiles
  #We exclude input probabilities of one, since quantiles for these are computed manually
  if (log.p) { LOGPROBS <- p } else { LOGPROBS <- log(p) }
  MAX.LOGP <- -Inf
  if (lower.tail) {
    for (i in 1:length(LOGPROBS)) {
      if (LOGPROBS[i] < 0)     { MAX.LOGP <- max(MAX.LOGP, LOGPROBS[i]) } } }
  if (!lower.tail) {
    for (i in 1:length(LOGPROBS)) {
      if (LOGPROBS[i] > - Inf) { MAX.LOGP <- max(MAX.LOGP, VGAM::log1mexp(-LOGPROBS[i])) } } }

  #Compute for special case where m = Inf
  if (m == Inf) {
    QUANTILES <- qnbinom(LOGPROBS, size = k, prob = prob, log.p = TRUE, lower.tail = lower.tail)
    return(QUANTILES) }

  #Compute quantiles using recursion
  if (!approx) {

    #Set log-probability matrix for t = 0
    t <- 0
    NEGOCC <- matrix(-Inf, nrow = 1, ncol = k+1)

    #Compute value for first row
    NEGOCC[1, 1] <- 0
    r <- 1
    while (r <= k) {
      NEGOCC[1, r+1] <- log(prob*(m-r+1)/m) + NEGOCC[1, r]
      r <- r+1 }
    CUML <- NEGOCC[1, k+1]

    #Iterate through t = 1,2,3,... until cumulative log-probability is large enough for quantiles
    while (CUML < MAX.LOGP) {

      #Increase value of t and add one row to NEGOCC
      t <- t+1
      NEGOCC <- rbind(NEGOCC, rep(-Inf, k+1))

      #Compute first value
      NEGOCC[t+1, 2] <- log(prob) + t*log(1-prob)

      #Compute values rows via recursion
      r <- 2
      while (r <= k) {
        LLL <- (0:t)*log(1-prob*(m-r+1)/m)
        TERMS <- LLL + NEGOCC[(t+1):1, r]
        NEGOCC[t+1, r+1] <- log(prob*(m-r+1)/m) + matrixStats::logSumExp(TERMS)
        r <- r+1 }

      #Update cumulative log-probability
      CUML <- matrixStats::logSumExp(c(CUML, NEGOCC[t+1, k+1])) }

    #Set maximum argument value
    max.x <- t

    #Generate cumulative log-probabilities
    CUMNEGOCC <- rep(-Inf, max.x+1)
    CUMNEGOCC[1] <- NEGOCC[1, k+1]
    if  (max.x > 0) {
    for (t in 1:max.x) {
      CUMNEGOCC[t+1] <- matrixStats::logSumExp(c(CUMNEGOCC[t], NEGOCC[t+1, k+1])) } }

    #Generate quantiles
    QUANTILES <- rep(0, length(LOGPROBS))
    for (i in 1:length(LOGPROBS)) {
      if (lower.tail)   { logprob <- LOGPROBS[i] } else { logprob <- VGAM::log1mexp(-LOGPROBS[i]) }
      if (logprob == 0) { QUANTILES[i] <- ifelse(k == 0, 0, ifelse(((k == 1)&(prob == 1)), 0, Inf)) }
      if (logprob < 0)  { QUANTILES[i] <- sum(CUMNEGOCC < logprob) } } }

  #Compute quantiles using approximation
  if (approx) {

    #Compute generalised harmonic numbers
    H1   <- sum(1/((m-k+1):m))
    H2   <- sum(1/((m-k+1):m)^2)

    #Compute moments and parameters
    MEAN <- max(0,(m/prob)*H1 - k)
    VAR  <- max(0,(m/prob)^2*H2 - (m/prob)*H1)
    SHAPE  <- (MEAN + 1/2)^2/VAR
    RATE   <- m*(MEAN + 1/2)/VAR

    #Set log-probabilities for t = 0
    t <- 0
    if (VAR == 0) {
      CUML <- 0
      CUMNEGOCC <- CUML }
    if (VAR > 0)  {
      CUML <- pgamma(1/m, shape = SHAPE, rate = RATE, log.p = TRUE)
      CUMNEGOCC <- CUML }

    #Iterate through t = 1,2,3,... until cumulative log-probability is large enough for quantiles
    while (CUML < MAX.LOGP) {

      #Approximation using discretised gamma distribution
      t <- t+1
      CUML <- pgamma((t+1)/m, shape = SHAPE, rate = RATE, log.p = TRUE)
      CUMNEGOCC <- c(CUMNEGOCC, CUML) }

    #Generate quantiles
    QUANTILES <- rep(0, length(LOGPROBS))
    for (i in 1:length(LOGPROBS)) {
      if (lower.tail)   { logprob <- LOGPROBS[i] } else { logprob <- VGAM::log1mexp(-LOGPROBS[i]) }
      if (logprob == 0) { QUANTILES[i] <- ifelse(k == 0, 0, ifelse(((k == 1)&(prob == 1)), 0, Inf)) }
      if (logprob < 0)  { QUANTILES[i] <- sum(CUMNEGOCC < logprob) } } }

  #Return output
  QUANTILES }
