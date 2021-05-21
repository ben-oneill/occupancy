#' Random generation from the negative occupancy distribution
#'
#' \code{rnegocc} returns random values from the distribution.
#'
#' This function generates random values from the negative occupancy distribution, which is the distribution for
#' the excess hitting time in the extended occupancy problem.
#'
#' @usage \code{rnegocc(n, space, occupancy, prob)}
#' @param n The number of random values to generate
#' @param space The space parameter for the negative occupancy distribution (number of bins)
#' @param occupancy The occupancy parameter for the negative occupancy distribution (number of occupied bins)
#' @param prob The probability parameter for the negative occupancy distribution (probability of ball occupying its bin)
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will be a
#' vector of random values of length \code{n}

rnegocc <- function(n, space = 1, occupancy = space, prob = 1) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(n))                        stop('Error: Argument n is not numeric')
  if (!is.numeric(space))                    stop('Error: Space parameter is not numeric')
  if (!is.numeric(occupancy))                stop('Error: Occupancy parameter is not numeric')
  if (!is.numeric(prob))                     stop('Error: Probability parameter is not numeric')

  #Check that parameters are atomic
  if (length(space)  != 1)                   stop('Error: Space parameter should be a single number')
  if (length(occupancy) != 1)                stop('Error: Occupancy parameter should be a single number')
  if (length(prob)  != 1)                    stop('Error: Probability parameter should be a single number')

  #Set parameters
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }
  k <- as.integer(occupancy)

  #Check that parameters are in allowable range
  if (length(n) != 1)                        stop('Error: Argument n should be a single positive integer')
  if (as.integer(n) != n)                    stop('Error: Argument n should be a positive integer')
  if (min(n) < 1)                            stop('Error: Argument n should be a positive integer')
  if (space != m)                            stop('Error: Space parameter is not an integer')
  if (m <= 0)                                stop('Error: Space parameter must be positive')
  if (occupancy != k)                        stop('Error: Occupancy parameter is not an integer')
  if (k > m)                                 stop('Error: Occupancy parameter is larger than space parameter')
  if (prob < 0)                              stop('Error: Probability parameter must be between zero and one')
  if (prob > 1)                              stop('Error: Probability parameter must be between zero and one')

  #Generate random values
  GEOM <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    if (m == Inf) { GEOM[, i] <- rgeom(n, prob = prob) }
    if (m < Inf)  { GEOM[, i] <- rgeom(n, prob = prob*(m-i+1)/m) } }
  OUT <- rowSums(GEOM)

  #Give output
  OUT }
