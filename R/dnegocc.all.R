#' Mass function of the negative occupancy distribution
#'
#' \code{dnegocc.all} returns array of probability or log-probability values up to a maximum argument.
#'
#' This function computes probabilities or log-probabilities from the mass function of the negative occupancy
#' distribution.  The computation method uses a recursive algorithm from the following paper:
#'
#' O'Neill, B. (2021) An examination of the negative-occupancy distribution and the coupon-collector distribution.
#'
#' @usage \code{dnegocc.all(max.x, space, occupancy, prob, approx = FALSE, log = FALSE)}
#' @param max.x A vector of numeric values to be used as arguments for the mass function
#' @param space The space parameter for the negative occupancy distribution (number of bins)
#' @param max.occupancy The maximum occupancy parameter for the negative occupancy distribution (number of occupied bins)
#' @param prob The probability parameter for the negative occupancy distribution (probability of ball occupying its bin)
#' @param approx A logical value specifying whether to use an approximation for the distribution
#' @param log A logical value specifying whether results should be returned as log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range and arguments are integers)
#' then the output will be a matrix of probabilities/log-probabilities

dnegocc.all <- function(max.x, space, max.occupancy, prob = 1, approx = FALSE, log = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(max.x))                    stop('Error: Argument max.x is not numeric')
  if (!is.numeric(space))                    stop('Error: Space parameter is not numeric')
  if (!is.numeric(max.occupancy))            stop('Error: Maximum occupancy parameter is not numeric')
  if (!is.numeric(prob))                     stop('Error: Probability parameter is not numeric')
  if (!is.logical(approx))                   stop('Error: approx option is not a logical value')
  if (!is.logical(log))                      stop('Error: log option is not a logical value')

  #Check that parameters are atomic
  if (length(max.x)  != 1)                   stop('Error: Argument max.x should be a single number')
  if (length(space)  != 1)                   stop('Error: Space parameter should be a single number')
  if (length(max.occupancy) != 1)            stop('Error: Maximum occupancy parameter should be a single number')
  if (length(prob)  != 1)                    stop('Error: Probability parameter should be a single number')
  if (length(approx) != 1)                   stop('Error: approx option should be a single logical value')
  if (length(log) != 1)                      stop('Error: log option should be a single logical value')

  #Set parameters
  mx <- as.integer(max.x)
  m  <- as.integer(space)
  k  <- as.integer(max.occupancy)

  #Check that parameters are in allowable range
  if (max.x != mx)                           stop('Error: Argument max.x is not an integer')
  if (mx < 0)                                stop('Error: Argument max.x must be non-negative')
  if (space != m)                            stop('Error: Space parameter is not an integer')
  if (m <= 0)                                stop('Error: Space parameter must be positive')
  if (max.occupancy != k)                    stop('Error: Maximum occupancy parameter is not an integer')
  if (k > m)                                 stop('Error: Maximum occupancy parameter is larger than space parameter')
  if (k < 0)                                 stop('Error: Maximum occupancy parameter is must be non-negative')
  if (prob < 0)                              stop('Error: Probability parameter must be between zero and one')
  if (prob > 1)                              stop('Error: Probability parameter must be between zero and one')

  #Create output vector
  NEGOCC <- matrix(-Inf, nrow = mx+1, ncol = k+1)
  rownames(NEGOCC) <- sprintf('t[%s]', 0:mx)
  colnames(NEGOCC) <- sprintf('k[%s]', 0:k)

  #Compute for trivial case where k = 0
  NEGOCC[1,1] <- 0
  if (k == 0) {
    if (log) { return(NEGOCC) } else { return(exp(NEGOCC)) } }

  #Compute for trivial case where prob = 0
  if (prob == 0) {
    if (log) { return(NEGOCC) } else { return(exp(NEGOCC)) } }

  #Compute for non-trivial cases where k > 0 and prob > 0
  #Compute log-probablities using recursion
  if (!approx) {

    #Compute first column of matrix
    if(prob == 1) {
      NEGOCC[ , 2] <- c(0, rep(-Inf, mx)) } else {
      NEGOCC[ , 2] <- log(prob) + (0:mx)*log(1-prob) }

    #Compute remaining rows via recursion
    if (k > 1) {
    for (r in 2:k) {
      LLL <- (0:mx)*log(1-prob*(m-r+1)/m)
      for (t in 0:mx) {
        TERMS <- LLL[1:(t+1)] + NEGOCC[(t+1):1,r]
        NEGOCC[t+1,r+1] <- log(prob*(m-r+1)/m) + matrixStats::logSumExp(TERMS) } } } }

  #Compute log-probabilities using approximation
  if (approx) {

    for (r in 1:k) {

      #Compute generalised harmonic numbers
      H1   <- sum(1/((m-r+1):m))
      H2   <- sum(1/((m-r+1):m)^2)

      #Compute moments
      MEAN <- max(0,(m/prob)*H1 - r)
      VAR  <- max(0,(m/prob)^2*H2 - (m/prob)*H1)

      #Approximation using discretised gamma distribution
      if (VAR == 0) {
        APPROX <- c(0, rep(-Inf, mx)) }
      if (VAR > 0)  {
        SHAPE  <- (MEAN + 1/2)^2/VAR
        RATE   <- m*(MEAN + 1/2)/VAR
        LGA    <- pgamma((0:(mx+1))/m, shape = SHAPE, rate = RATE, log.p = TRUE)
        LOWER  <- LGA[1:(mx+1)]
        UPPER  <- LGA[2:(mx+2)]
        APPROX <- UPPER + VGAM::log1mexp(UPPER-LOWER)
        NEGOCC[, r+1] <- APPROX } } }

  #Return output
  if (log) { NEGOCC } else { exp(NEGOCC) } }
