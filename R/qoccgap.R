#' @rdname doccgap
qoccgap <- function(p, size, space = NULL, occupancy = size, prob = NULL, scale = NULL, log.p = FALSE, lower.tail = TRUE) {

  #Check scale parameter
  if (!is.null(scale)) {
    if (!is.numeric(scale))                    stop('Error: Scale parameter is not numeric')
    if (length(scale) != 1)                    stop('Error: Scale parameter should be a single number')
    if (scale < 0)                             stop('Error: Scale parameter must be non-negative') }

  #Check space parameter
  if (!is.null(space)) {
    if (!is.numeric(space))                    stop('Error: Space parameter is not numeric')
    if (length(space) != 1)                    stop('Error: Space parameter should be a single number')
    m <- as.integer(space)
    if (space != m)                            stop('Error: Size parameter should be a single number')
    if (m < 0)                                 stop('Error: Space parameter must be non-negative') }

  #Check probability parameter
  if (!is.null(prob)) {
    if (!is.numeric(prob))                     stop('Error: Probability parameter is not numeric')
    if (length(prob) != 1)                     stop('Error: Probability parameter should be a single number')
    if (prob < 0)                              stop('Error: Probability parameter must be between zero and one')
    if (prob > 1)                              stop('Error: Probability parameter must be between zero and one') }

  #Check parameterisation
  if (!is.null(scale)) {
    if ((!is.null(space))&(is.null(prob)))     stop('Error: Specify scale parameter or space and probability, but not both')
    if ((is.null(space))&(!is.null(prob)))     stop('Error: Specify scale parameter or space and probability, but not both')
    if ((!is.null(space))&(!is.null(prob))) {
      ERR <- abs(scale - m*(1-prob)/prob)
      if (ERR <= 1e-6) {
        warning('Specify scale parameter or space and probability, but not both') } else {
        stop('Error: Specify scale parameter or space and probability, but not both') } } }
  if (is.null(scale)) {
    if ((is.null(space))|(is.null(prob)))      stop('Error: You must either specify scale parameter or space and probability')
    scale <- m*(1-prob)/prob }

  #Check that argument and parameters are appropriate type
  if (!is.numeric(p))                        stop('Error: Argument p is not numeric')
  if (!is.numeric(size))                     stop('Error: Size parameter is not numeric')
  if (!is.numeric(occupancy))                stop('Error: Occupancy parameter is not numeric')
  if (!is.logical(log.p))                    stop('Error: log.p option is not a logical value')
  if (!is.logical(lower.tail))               stop('Error: lower.tail option is not a logical value')

  #Check that parameters are atomic
  if (length(size)  != 1)                    stop('Error: Size parameter should be a single number')
  if (length(occupancy) != 1)                stop('Error: Occupancy parameter should be a single number')
  if (length(log.p) != 1)                    stop('Error: log.p option should be a single logical value')
  if (length(lower.tail) != 1)               stop('Error: lower.tail option should be a single logical value')

  #Check that argument values are in allowable range
  if (!log.p) {
    if (min(p) < 0)                          stop('Error: Probability values in p must be between zero and one')
    if (max(p) > 1)                          stop('Error: Probability values in p must be between zero and one') }
  if (log.p) {
    if (max(p) > 0)                          stop('Error: Log-probability values in p must be less than or equal to zero') }

  #Set parameters
  n <- as.integer(size)
  k <- as.integer(occupancy)

  #Check that parameters are in allowable range
  if (size != n)                             stop('Error: Size parameter is not an integer')
  if (n <= 0)                                stop('Error: Size parameter must be positive')
  if (occupancy != k)                        stop('Error: Occupancy parameter is not an integer')
  if (k < 0)                                 stop('Error: Occupancy parameter is must be non-negative')
  if (k > n)                                 stop('Error: Occupancy parameter is larger than size parameter')
  if (!is.null(space)) {
    if (k > m)                               stop('Error: Occupancy parameter is larger than space parameter') }

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

  #Compute log-probablities using recursion
    #Set log-Stirling matrix and generate first row
    LOGSTIRLING <- matrix(-Inf, nrow = n+1, ncol = k+1)
    LOGSTIRLING[1,1] <- 0

    #Generate subsequent rows
    if (k > 0) {
     for (nn in 1:n) {
      for (kk in 1:min(k,nn)) {
        T1 <- log(kk) + LOGSTIRLING[nn, kk+1]
        T2 <- LOGSTIRLING[nn, kk]
        LOGSTIRLING[nn+1, kk+1] <- matrixStats::logSumExp(c(T1, T2)) } } }

    #Generate the cumulative log-probabilities for the occupancy-gap distribution
    LOGS    <- rep(-Inf, n-k+1)
    CUMLOGS <- rep(-Inf, n-k+1)
    for (i in k:n) {
      LOGS[i-k+1]    <- lchoose(n,i) + (n-i)*log(scale) + LOGSTIRLING[i+1, k+1] }
    LOGS <- LOGS - matrixStats::logSumExp(LOGS)
    for (i in k:n) {
      CUMLOGS[i-k+1] <- matrixStats::logSumExp(LOGS[1:(i-k+1)]) }

    #Generate quantiles
    QUANTILES <- rep(0, length(LOGPROBS))
    for (i in 1:length(LOGPROBS)) {
      if (lower.tail)   { logprob <- LOGPROBS[i] } else { logprob <- VGAM::log1mexp(-LOGPROBS[i]) }
      QUANTILES[i] <- sum(CUMLOGS < logprob) }

  #Return output
  QUANTILES }
