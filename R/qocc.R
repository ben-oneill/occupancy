#' @rdname docc
qocc <- function(p, size, space, prob = 1, approx = FALSE, log.p = FALSE, lower.tail = TRUE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(p))                                  stop('Error: Argument p is not numeric')
  if (!is.numeric(size))                               stop('Error: Size parameter is not numeric')
  if (!is.numeric(space))                              stop('Error: Space parameter is not numeric')
  if (!is.numeric(prob))                               stop('Error: Probability parameter is not numeric')
  if (!is.logical(approx))                             stop('Error: approx option is not a logical value')
  if (!is.logical(log.p))                              stop('Error: log.p option is not a logical value')
  if (!is.logical(lower.tail))                         stop('Error: lower.tail option is not a logical value')

  #Check that parameters are atomic
  if (length(size)  != 1)                              stop('Error: Size parameter should be a single number')
  if (length(space) != 1)                              stop('Error: Space parameter should be a single number')
  if (length(prob)  != 1)                              stop('Error: Probability parameter should be a single number')
  if (length(approx)   != 1)                           stop('Error: approx option should be a single logical value')
  if (length(log.p) != 1)                              stop('Error: log.p option should be a single logical value')
  if (length(lower.tail) != 1)                         stop('Error: lower.tail option should be a single logical value')

  #Set parameters
  n   <- as.integer(size)
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }
  MAX <- min(n,m)

  #Check that parameters are in allowable range
  if (size != n)                                       stop('Error: Size parameter is not an integer')
  if (n < 0)                                           stop('Error: Size parameter must be non-negative')
  if (space != m)                                      stop('Error: Space parameter is not an integer')
  if (m <= 0)                                          stop('Error: Space parameter must be positive')
  if ((prob < 0)|(prob > 1))                           stop('Error: Probability parameter is not between zero and one')

  #Check that argument values are in allowable range
  if (!log.p) {
    if (min(p) < 0)                                    stop('Error: Probability values in p must be between zero and one')
    if (max(p) > 1)                                    stop('Error: Probability values in p must be between zero and one') }
  if (log.p) {
    if (max(p) > 0)                                    stop('Error: Log-probability values in p must be less than or equal to zero') }

  #Compute log-probabilities for input p
  #Adjust for floating point error in computation log(exp(...))
  MAX.FP.ERROR <- .Machine$double.eps
  if (lower.tail) {
    if (log.p) { LOGP <- p } else { LOGP <- log(p) - MAX.FP.ERROR } }
  if (!lower.tail) {
    if (log.p) { LOGP <- VGAM::log1mexp(-p) } else { LOGP <- log(1-p) - MAX.FP.ERROR } }

  #Compute for trivial case where n = 0 or prob = 0
  if ((n == 0)|(prob == 0)) {
    QUANTILES <- rep(0, length(p))
    return(QUANTILES) }

  #Compute for special case where m = Inf
  if (m == Inf) {
    QUANTILES <- qbinom(LOGPROBS, size = n, prob = prob, log.p = TRUE, lower.tail = lower.tail)
    return(QUANTILES) }

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

  #Generate quantiles
  QUANTILES <- rep(0, length(LOGPROBS))
  for (i in 1:length(LOGPROBS)) {
    logp <- LOGP[i]
    if (logp == 0) { QUANTILES[i] <- MAX } else { QUANTILES[i] <- sum(CUMLOGS < logp) } }

  #Return output
  QUANTILES }
