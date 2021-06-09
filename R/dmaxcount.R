#' @rdname dmaxcount.all
dmaxcount <- function(x, size, space, prob = 1, log = FALSE) {

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
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }

  #Check that parameters are in allowable range
  if (size != n)                stop('Error: Size parameter is not an integer')
  if (n < 0)                    stop('Error: Size parameter should be nonnegative')
  if (space != m)               stop('Error: Space parameter is not an integer')
  if (m <= 0)                   stop('Error: Space parameter should be positive')
  if ((prob < 0)|(prob > 1))    stop('Error: Probability parameter is not between zero and one')

  ################################################################################################################
  ######  Compute the cumulative log-probabilities via iterative method in Bonetti, Corillo and Ogay (2019)  #####
  ################################################################################################################

  #Deal with trivial case where n = 0
  if (n == 0) {
    OUT <- rep(-Inf, length(x))
    IND <- (x == 0)
    OUT[IND] <- 0
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with trival case where m = Inf
  if (m == Inf) {

    #Compute probabilities of max-count = 0 or 1
    P0 <- n*log(1-prob)
    P1 <- VGAM::log1mexp(-P0)

    #Compute and return output vector
    OUT <- rep(-Inf, length(x))
    for (i in 1:length(x)) {
      xx <- x[i]
      if (xx == 0) { OUT[i] <- P0 }
      if (xx == 1) { OUT[i] <- P1 } }
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Deal with non-trivial case where m < Inf
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
    LOGPROBS   <- LLL[, n+1, m]
  } else {
    LOGBINDIST <- dbinom(x = 0:n, size = n, prob = prob, log = TRUE)
    LOGPROBS   <- rep(-Inf, MAX+1)
    for (xx in 0:MAX) {
      LOGPROBS[xx+1] <- matrixStats::logSumExp(LLL[xx+1, , m] + LOGBINDIST) } }

  #Compute the function output (the log-probabilities over the argument values)
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
