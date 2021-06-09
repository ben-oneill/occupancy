#' @rdname dnegocc.all
dnegocc <- function(x, space, occupancy, prob = 1, approx = FALSE, log = FALSE) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(x))                        stop('Error: Argument x is not numeric')
  if (!is.numeric(space))                    stop('Error: Space parameter is not numeric')
  if (!is.numeric(occupancy))                stop('Error: Occupancy parameter is not numeric')
  if (!is.numeric(prob))                     stop('Error: Probability parameter is not numeric')
  if (!is.logical(approx))                   stop('Error: approx option is not a logical value')
  if (!is.logical(log))                      stop('Error: log option is not a logical value')

  #Check that parameters are atomic
  if (length(space)  != 1)                   stop('Error: Space parameter should be a single number')
  if (length(occupancy) != 1)                stop('Error: Occupancy parameter should be a single number')
  if (length(prob)  != 1)                    stop('Error: Probability parameter should be a single number')
  if (length(approx) != 1)                   stop('Error: approx option should be a single logical value')
  if (length(log) != 1)                      stop('Error: log option should be a single logical value')

  #Set parameters
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }
  k <- as.integer(occupancy)

  #Check that parameters are in allowable range
  if (space != m)                            stop('Error: Space parameter is not an integer')
  if (m <= 0)                                stop('Error: Space parameter must be positive')
  if (occupancy != k)                        stop('Error: Occupancy parameter is not an integer')
  if (k < 0)                                 stop('Error: Occupancy parameter must be non-negative')
  if (k > m)                                 stop('Error: Occupancy parameter is larger than space parameter')
  if (prob < 0)                              stop('Error: Probability parameter must be between zero and one')
  if (prob > 1)                              stop('Error: Probability parameter must be between zero and one')

  #Create output vector
  max.x  <- floor(max(x))
  NEGOCC <- rep(-Inf, length(x))

  #Compute for trivial case where k = 0
  if (k == 0) {
    for (i in 1:length(x)) {
      xx <- x[i]
      if (xx == 0) { NEGOCC[i] <- 0 } }
    if (log) { return(NEGOCC) } else { return(exp(NEGOCC)) } }

  #Compute for trivial case where prob = 0
  if (prob == 0) {
  if (k > 0) {
    for (i in 1:length(x)) {
      xx <- x[i]
      if (xx == Inf) { NEGOCC[i] <- 0 } }
    if (log) { return(NEGOCC) } else { return(exp(NEGOCC)) } } }

  #Compute for special case where m = Inf
  if (m == Inf) {
    NEGOCC <- dnbinom(x, size = k, prob = prob, log = TRUE)
    if (log) { return(NEGOCC) } else { return(exp(NEGOCC)) } }

  #Compute for non-trivial cases where k > 0 and prob > 0
  #Compute log-probablities using recursion
  if (!approx) {

    #Create base vector for recursion
    if(prob == 1) {
      LOGS <- c(0, rep(-Inf, max.x)) } else {
      LOGS <- log(prob) + (0:max.x)*log(1-prob) }

    #Update via recursion
    r <- 2
    while (r <= k) {
      NEWLOGS <- rep(-Inf, max.x+1)
      LLL <- (0:max.x)*log(1-prob*(m-r+1)/m)
      for (t in 0:max.x) {
        TERMS <- LLL[1:(t+1)] + LOGS[(t+1):1]
        NEWLOGS[t+1] <- log(prob*(m-r+1)/m) + matrixStats::logSumExp(TERMS) }
      LOGS <- NEWLOGS
      r    <- r+1 }

    #Generate output vector
    for (i in 1:length(x)) {
      xx <- x[i]
      if ((as.integer(xx) == xx)&(xx >= 0)) {
        NEGOCC[i] <- LOGS[xx+1] } } }

  #Compute log-probabilities using approximation
  if (approx) {

    #Compute generalised harmonic numbers
    H1   <- sum(1/((m-k+1):m))
    H2   <- sum(1/((m-k+1):m)^2)

    #Compute moments
    MEAN <- max(0,(m/prob)*H1 - k)
    VAR  <- max(0,(m/prob)^2*H2 - (m/prob)*H1)

    #Approximation using discretised gamma distribution
    if (VAR == 0) {
      APPROX <- c(0, rep(-Inf, max.x)) }
    if (VAR > 0)  {
      SHAPE  <- (MEAN + 1/2)^2/VAR
      RATE   <- m*(MEAN + 1/2)/VAR
      LGA    <- pgamma((0:(max.x+1))/m, shape = SHAPE, rate = RATE, log.p = TRUE)
      LOWER  <- LGA[1:(max.x+1)]
      UPPER  <- LGA[2:(max.x+2)]
      APPROX <- UPPER + VGAM::log1mexp(UPPER-LOWER) }

    #Generate output vector
    for (i in 1:length(x)) {
      xx <- x[i]
      if ((as.integer(xx) == xx)&(xx >= 0)) {
        NEGOCC[i] <- APPROX[xx+1] } } }

  #Return output
  if (log) { NEGOCC } else { exp(NEGOCC) } }
