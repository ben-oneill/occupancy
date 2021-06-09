#' @rdname doccgap.all
doccgap <- function(x, size, space = NULL, occupancy = size, prob = NULL, scale = NULL, log = FALSE) {

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
  if (!is.numeric(x))                          stop('Error: Argument x is not numeric')
  if (!is.numeric(size))                       stop('Error: Size parameter is not numeric')
  if (!is.numeric(occupancy))                  stop('Error: Occupancy parameter is not numeric')
  if (!is.logical(log))                        stop('Error: log option is not a logical value')

  #Check that parameters are atomic
  if (length(size)  != 1)                      stop('Error: Size parameter should be a single number')
  if (length(occupancy) != 1)                  stop('Error: Occupancy parameter should be a single number')
  if (length(log) != 1)                        stop('Error: log option should be a single logical value')

  #Set parameters
  n <- as.integer(size)
  k <- as.integer(occupancy)

  #Check that parameters are in allowable range
  if (size != n)                               stop('Error: Size parameter is not an integer')
  if (n < 0)                                   stop('Error: Size parameter must be non-negative')
  if (occupancy != k)                          stop('Error: Occupancy parameter is not an integer')
  if (k < 0)                                   stop('Error: Occupancy parameter is must be non-negative')
  if (k > n)                                   stop('Error: Occupancy parameter is larger than size parameter')
  if (!is.null(space)) {
    if (k > m)                                 stop('Error: Occupancy parameter is larger than space parameter') }

  #Deal with trivial case where n = 0
  if (n == 0) {
    OUT <- rep(-Inf, length(x))
    IND <- (x == 0)
    OUT[IND] <- 0
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Create output vector
  OCCGAP <- rep(-Inf, length(x))

  #Compute for trivial case where scale = 0
  if (scale == 0) {
    for (i in 1:length(x)) {
      xx <- x[i]
      if (xx == n-k) { OCCGAP[i] <- 0 } }
    if (log) { return(OCCGAP) } else { return(exp(OCCGAP)) } }

  #Compute for trivial case where scale = Inf
  if (scale == Inf) {
    for (i in 1:length(x)) {
      xx <- x[i]
      if (xx == 0) { OCCGAP[i] <- 0 } }
    if (log) { return(OCCGAP) } else { return(exp(OCCGAP)) } }

  #Compute for non-trivial cases where 0 < scale < Inf
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

  #Generate the log-probabilities for the occupancy-gap distribution
  LOGS <- rep(-Inf, n-k+1)
  for (i in k:n) {
    LOGS[i-k+1] <- lchoose(n,i) + (n-i)*log(scale) + LOGSTIRLING[i+1, k+1] }
  LOGS <- LOGS - matrixStats::logSumExp(LOGS)

  #Generate output vector
  for (i in 1:length(x)) {
    xx <- x[i]
    if (xx %in% 0:(n-k)) {
      OCCGAP[i] <- LOGS[xx+1] } }

  #Return output
  if (log) { OCCGAP } else { exp(OCCGAP) } }
