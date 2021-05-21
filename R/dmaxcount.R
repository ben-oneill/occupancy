#' Probabilty mass function of the maximum-count distribution
#'
#' \code{dmaxcount} returns the probability or log-probability values for the arguments.
#'
#' This function computes probabilities or log-probabilities from the probability mass function of the maximum-count
#' distribution, which is the distribution for the maximum of the counts for the number of balls in a bin in the extended
#' occupancy problem.  Details of the algorithm in the classical case can be found in the papers below.  The extension
#' to include the probability parameter is done using the binomial mixture representation of the extended occupancy problem.
#'
#' Bonetti, M. and Corillo, P. (2019) Computing the exact distributions of some functions of the ordered multinomial
#' counts: maximum, minimum, range and sums of order statistics.
#'
#' Rappeport, M,A. (1968) Algorithms and computational procedures for the application of order statistics to queuing
#' problems. PhD thesis, New York University.
#'
#' @usage \code{dmaxcount(x, size, space, prob, log = FALSE)}
#' @param x A vector of numeric values to be used as arguments for the probability mass function
#' @param size The size parameter for the maximum-count distribution (number of balls)
#' @param space The space parameter for the maximum-count distribution (number of bins)
#' @param prob The probability parameter for the occupancy distribution (probability of ball occupying its bin)
#' @param log A logical value specifying whether results should be returned as log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will be a
#' vector of probabilities/log-probabilities corresponding to the vector argument x

dmaxcount <- function(x, size = 1, space = 1, prob = 1, log = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(x))           stop('Error: Argument x is not numeric')
  if (!is.numeric(size))        stop('Error: Size parameter is not numeric')
  if (!is.numeric(space))       stop('Error: Space parameter is not numeric')
  if (!is.numeric(prob))        stop('Error: Probability parameter is not numeric')
  if (!is.logical(log))         stop('Error: log option is not a logical value')

  #Check that parameters are atomic
  if (length(size)  != 1)       stop('Error: Size parameter should be a single number')
  if (length(space) != 1)       stop('Error: Space parameter should be a single number')
  if (length(prob)  != 1)       stop('Error: Probability parameter should be a single number')
  if (length(log) != 1)         stop('Error: log.p option should be a single logical value')

  #Set parameters
  n <- as.integer(size)
  m <- as.integer(space)

  #Check that parameters are in allowable range
  if (size != n)                stop('Error: Size parameter is not an integer')
  if (n <= 0)                   stop('Error: Size parameter is negative')
  if (space != m)               stop('Error: Space parameter is not an integer')
  if (m <= 0)                   stop('Error: Space parameter is negative')
  if ((prob < 0)|(prob > 1))    stop('Error: Probability parameter is not between zero and one')

  ################################################################################################################
  #########  Compute the cumulative log-probabilities via iterative method in Bonetti and Corillo (2019)  ########
  ################################################################################################################

  #Create matrix of log-probabilities
  MAX <- min(floor(max(x)), n)
  LLL <- array(-Inf, dim = c(MAX+1, n+1, m),
               dimnames = list(sprintf('t[%s]', 0:MAX), sprintf('n[%s]', 0:n), sprintf('m[%s]', 1:m)))

  #Set trivial log-probabilities in the case where nn <= xx (i.e., no more balls than xx)
  for (xx in 0:MAX) { LLL[xx+1, 1:(xx+1), ] <- 0 }

  #Compute remaining non-trivial log-probabilities
  if (m > 1) {

    #Compute the base log-probabilities for maxcount t = 1
    for (nn in 1:n) {
    for (mm in 2:m) {
      if (mm >= nn) { LLL[2, nn+1, mm] <- lchoose(mm, nn) + lfactorial(nn) - nn*log(mm) } } }

    #Iteratively compute the log-probabilities for maxcount t > 1
    if (MAX > 1) {
    for (xx in 2:MAX) {
      for (nn in 1:n)  {
      for (mm in 2:m)  {
      if (nn > xx) {

        #Generate weighting terms
        ITER  <- FALSE
        LOWER <- max(0, nn-xx*mm+mm)
        UPPER <- floor(nn/xx)
        if (UPPER >= LOWER) {
          QQ   <- LOWER:UPPER
          ITER <- TRUE }

        #Compute terms for iteration
        if (ITER) {

          #Set terms for recursion
          LOGA  <- rep(-Inf, length(QQ))
          LOGP  <- rep(-Inf, length(QQ))
          for (i in 1:length(QQ)) {
            qq <- QQ[i]
            if ((nn >= xx*qq)&(mm > qq)) {
              LOGA[i] <- lfactorial(nn) + lfactorial(mm) - nn*log(mm) - qq*lfactorial(xx) - lfactorial(qq) +
                         (nn-xx*qq)*log(mm-qq) - lfactorial(mm-qq) - lfactorial(nn-xx*qq)
              LOGP[i] <- LLL[xx, nn-xx*qq+1, mm-qq] }
            if ((nn == xx*qq)&(mm == qq)) {
              LOGA[i] <- lfactorial(nn) + lfactorial(mm) - nn*log(mm) - qq*lfactorial(xx) - lfactorial(qq)
              LOGP[i] <- 0 } }

        #Compute the new log-probability
          LLL[xx+1, nn+1, mm] <- min(0, max(matrixStats::logSumExp(LOGA + LOGP), LLL[xx, nn+1, mm])) } } } } } } }

  #Compute the log-probabilities for the maximum count for the extended occupancy problem
  if (prob == 1) {
    LOGPROBS   <- LLL[, n+1, m] } else {
    LOGBINDIST <- dbinom(x = 0:n, size = n, prob = prob, log = TRUE)
    LOGPROBS   <- rep(-Inf, MAX+1)
    for (xx in 0:MAX) {
      LOGPROBS[xx+1] <- matrixStats::logSumExp(LLL[xx+1, , m] + LOGBINDIST) } }

  #Compute the function output (the log-CDF over the argument values)
  OUT <- rep(-Inf, length(x))
  for (i in 1:length(x)) {
    xx <- x[i]
    if (xx == 0) {
      OUT[i] <- LOGPROBS[1] }
    if (xx %in% 1:n) {
      if (LOGPROBS[xx+1] > LOGPROBS[xx]) {
        OUT[i] <- LOGPROBS[xx+1] + VGAM::log1mexp(LOGPROBS[xx+1] - LOGPROBS[xx]) } } }

  #Return the output
  if (log) { OUT } else { exp(OUT) } }
