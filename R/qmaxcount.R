#' @rdname dmaxcount
qmaxcount <- function(p, size, space, prob = 1, log.p = FALSE, lower.tail = TRUE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(p))                       stop('Error: Argument p is not numeric')
  if (!is.numeric(size))                    stop('Error: Size parameter is not numeric')
  if (!is.numeric(space))                   stop('Error: Space parameter is not numeric')
  if (!is.numeric(prob))                    stop('Error: Probability parameter is not numeric')
  if (!is.logical(log.p))                   stop('Error: log.p option is not a logical value')
  if (!is.logical(lower.tail))              stop('Error: lower.tail option is not a logical value')

  #Check that parameters are atomic
  if (length(size)  != 1)                   stop('Error: Size parameter should be a single number')
  if (length(space) != 1)                   stop('Error: Space parameter should be a single number')
  if (length(prob)  != 1)                   stop('Error: Probability parameter should be a single number')
  if (length(log.p) != 1)                   stop('Error: log.p option should be a single logical value')
  if (length(lower.tail) != 1)              stop('Error: lower.tail option should be a single logical value')

  #Set parameters
  n <- as.integer(size)
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }

  #Check that parameters are in allowable range
  if (size != n)                            stop('Error: Size parameter is not an integer')
  if (n < 0)                                stop('Error: Size parameter should be non-negative')
  if (space != m)                           stop('Error: Space parameter is not an integer')
  if (m <= 0)                               stop('Error: Space parameter should be positive')
  if ((prob < 0)|(prob > 1))                stop('Error: Probability parameter is not between zero and one')

  #Check that argument values are in allowable range
  if (!log.p) {
    if (min(p) < 0)                         stop('Error: Probability values in p must be between zero and one')
    if (max(p) > 1)                         stop('Error: Probability values in p must be between zero and one') }
  if (log.p) {
    if (max(p) > 0)                         stop('Error: Log-probability values in p must be less than or equal to zero') }

  ################################################################################################################
  #########  Compute the cumulative log-probabilities via iterative method in Bonetti and Corillo (2019)  ########
  ################################################################################################################

  #Deal with the trivial case where n = 0
  if (n == 0) {
    QUANTILE <- rep(0, length(p))
    return(QUANTILE) }

  #Create matrix of log-probabilities
  LLL  <- array(-Inf, dim = c(n+1, n+1, m),
                dimnames = list(sprintf('t[%s]', 0:n), sprintf('n[%s]', 0:n), sprintf('m[%s]', 1:m)))

  #Set trivial log-probabilities in the case where nn <= xx (i.e., no more balls than xx)
  for (xx in 0:n) { LLL[xx+1, 1:(xx+1), ] <- 0 }

  #Compute remaining non-trivial log-probabilities
  if (m > 1) {

    #Compute the base log-probabilities for maxcount t = 1
    for (nn in 1:n) {
      for (mm in 2:m) {
        if (mm >= nn) { LLL[2, nn+1, mm] <- lchoose(mm, nn) + lfactorial(nn) - nn*log(mm) } } }

    #Iteratively compute the log-probabilities for maxcount t > 1
    if (n > 1) {
      for (xx in 2:n) {
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
      LOGPROBS   <- rep(-Inf, n+1)
      for (xx in 0:n) {
        LOGPROBS[xx+1] <- matrixStats::logSumExp(LLL[xx+1, , m] + LOGBINDIST) } }

  #Obtain log-probabilities for quantiles
  if (lower.tail) {
    if (log.p) { LOGP <- p } else { LOGP <- log(p) } }
  if (!lower.tail) {
    if (log.p) { LOGP <- VGAM::log1mexp(-p) } else { LOGP <- log(1-p) } }

  #Compute quantiles
  RR <- length(p)
  QUANTILE <- rep(NA, RR)
  for (i in 1:RR) {
    Q <- 0
    L <- LOGPROBS[Q+1]
    while (LOGP[i] > L) {
      Q <- Q+1
      L <- LOGPROBS[Q+1] }
    QUANTILE[i] <- Q }

  #Return the output
  QUANTILE }
