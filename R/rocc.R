#' @rdname docc
rocc <- function(n, size, space, prob = 1) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(n))                       stop('Error: Argument K is not numeric')
  if (!is.numeric(size))                    stop('Error: Size parameter is not numeric')
  if (!is.numeric(space))                   stop('Error: Space parameter is not numeric')
  if (!is.numeric(prob))                    stop('Error: Probability parameter is not numeric')

  #Check that parameters are atomic
  if (length(size)  != 1)                   stop('Error: Size parameter should be a single number')
  if (length(space) != 1)                   stop('Error: Space parameter should be a single number')
  if (length(prob)  != 1)                   stop('Error: Probability parameter should be a single number')

  #Set parameters
  obs <- n
  n   <- as.integer(size)
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }
  MAX <- min(n,m)

  #Check that parameters are in allowable range
  if (size != n)                            stop('Error: Size parameter is not an integer')
  if (n < 0)                                stop('Error: Size parameter should be non-negative')
  if (space != m)                           stop('Error: Space parameter is not an integer')
  if (m <= 0)                               stop('Error: Space parameter is negative')
  if ((prob < 0)|(prob > 1))                stop('Error: Probability parameter is not between zero and one')

  #Generate random occupancy values
  ROCC   <- rep(0, obs)
  for (i in 1:obs) {
    if (n == 0) {
      ROCC[i] <- 0 }
    if (n > 0) {
      SAMPLE <- sample.int(n = m, size = n, replace = TRUE)
      IND    <- (runif(n) <= prob)
      SAMPLE[!IND] <- NA
      ROCC[i] <- length(table(SAMPLE)) } }

  #Output the vector of random values
  ROCC }
