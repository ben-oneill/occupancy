#' Probabilty mass function of the maximum-count distribution
#'
#' \code{dmaxcount.all} returns a matrix of probability or log-probability values up to a maximum arguments.
#'
#' This function computes probabilities or log-probabilities from the probability mass function of the maximum-count
#' distribution, which is the distribution for the maximum of the counts for the number of balls in a bin in the extended
#' occupancy problem.  Details of the algorithm in the classical case can be found in the papers below.  The extension
#' to include the probability parameter is done using the binomial mixture representation of the extended occupancy problem.
#'
#' Bonetti, M., Corillo, P. and Ogay, A. (2019) Computing the exact distributions of some functions of the ordered multinomial
#' counts: maximum, minimum, range and sums of order statistics.
#'
#' Rappeport, M,A. (1968) Algorithms and computational procedures for the application of order statistics to queuing
#' problems. PhD thesis, New York University.
#'
#' @param max.x A vector of numeric values to be used as arguments for the probability mass function
#' @param max.size The maximum size parameter for the maximum-count distribution (number of balls)
#' @param space The space parameter for the maximum-count distribution (number of bins)
#' @param prob The probability parameter for the occupancy distribution (probability of ball occupying its bin)
#' @param log A logical value specifying whether results should be returned as log-probabilities
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will be a
#' vector of probabilities/log-probabilities up to the maximum argument values

dmaxcount.all <- function(max.x, max.size, space, prob = 1, log = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(max.x))                   stop('Error: Argument max.x is not numeric')
  if (!is.numeric(max.size))                stop('Error: Maximum size parameter is not numeric')
  if (!is.numeric(space))                   stop('Error: Space parameter is not numeric')
  if (!is.numeric(prob))                    stop('Error: Probability parameter is not numeric')
  if (!is.logical(log))                     stop('Error: log option is not a logical value')

  #Check that parameters are atomic
  if (length(max.x)  != 1)                  stop('Error: Argument max.x should be a single number')
  if (length(max.size) != 1)                stop('Error: Maximum size parameter should be a single number')
  if (length(space) != 1)                   stop('Error: Apace parameter should be a single number')
  if (length(prob) != 1)                    stop('Error: Probability parameter should be a single number')
  if (length(log) != 1)                     stop('Error: log.p option should be a single logical value')

  #Set parameters
  n <- as.integer(max.size)
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }

  #Check that parameters are in allowable range
  if (size != n)                            stop('Error: Size parameter is not an integer')
  if (n < 0)                                stop('Error: Size parameter should be nonnegative')
  if (space != m)                           stop('Error: Space parameter is not an integer')
  if (m <= 0)                               stop('Error: Space parameter should be positive')
  if ((prob < 0)|(prob > 1))                stop('Error: Probability parameter is not between zero and one')

  ################################################################################################################
  ######  Compute the cumulative log-probabilities via iterative method in Bonetti, Corillo and Ogay (2019)  #####
  ################################################################################################################

  #Deal with trivial case where n = 0
  if (n == 0) {
    MAXCOUNT <- matrix(-Inf, nrow = MAX+1, ncol = 1)
    rownames(MAXCOUNT) <- sprintf('t[%s]', 0:MAX)
    colnames(MAXCOUNT) <- 'n[0]'
    MAXCOUNT[1, 1] <- 0
    if (log) { return(MAXCOUNT) } else { return(exp(MAXCOUNT)) } }

  #Deal with trival case where m = Inf
  if (m == Inf) {

    #Create output matrix
    MAX <- max.x
    MAXCOUNT <- matrix(-Inf, nrow = MAX+1, ncol = n+1)
    rownames(MAXCOUNT) <- sprintf('t[%s]', 0:MAX)
    colnames(MAXCOUNT) <- sprintf('n[%s]', 0:n)

    #Compute and return output matrix
    MAXCOUNT[1, 1] <- 0
    for (nn in 1:n) {
      P0 <- nn*log(1-prob)
      P1 <- VGAM::log1mexp(-P0)
      MAXCOUNT[1, nn+1] <- P0
      if (MAX > 0) { MAXCOUNT[2, nn+1] <- P1 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with non-trivial case where m < Inf
  #Create matrix of log-probabilities
  MAX <- max.x
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
    LOGPROBS <- LLL }
  if (prob < 1) {
    LOGPROBS <- array(-Inf, dim = c(MAX+1, n+1, m),
                      dimnames = list(sprintf('t[%s]', 0:MAX), sprintf('n[%s]', 0:n), sprintf('m[%s]', 1:m)))
    for (nn in 0:n) {
      LOGBINDIST <- dbinom(x = 0:nn, size = nn, prob = prob, log = TRUE)
      for (xx in 0:MAX) {
      for (mm in 1:m) {
        LOGPROBS[xx+1, nn+1, m] <- matrixStats::logSumExp(LLL[xx+1, 1:(nn+1), m] + LOGBINDIST) } } } }

  #Compute the log-probabilities for the mass function
  MAXCOUNT <- matrix(-Inf, nrow = MAX+1, ncol = n+1)
  rownames(MAXCOUNT) <- sprintf('t[%s]', 0:MAX)
  colnames(MAXCOUNT) <- sprintf('n[%s]', 0:n)
  MAXCOUNT[1, 1] <- 0
  for (nn in 1:n) {
  for (xx in 0:MAX) {
    if (xx == 0) {
      MAXCOUNT[xx+1, nn+1] <- LOGPROBS[xx+1, nn+1, m] }
    if (xx > 0) {
      L1 <- LOGPROBS[xx+1, nn+1, m]
      L0 <- LOGPROBS[xx, nn+1, m]
      if (L1 > L0) {
        MAXCOUNT[xx+1, nn+1] <- L1 + VGAM::log1mexp(L1 - L0) } } } }

  #Return the output
  if (log) { MAXCOUNT } else { exp(MAXCOUNT) } }
