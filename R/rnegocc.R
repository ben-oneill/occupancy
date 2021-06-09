#' @rdname dnegocc.all
rnegocc <- function(n, space, occupancy, prob = 1) {

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
  if (k < 0)                                 stop('Error: Occupancy parameter must be non-negative')
  if (k > m)                                 stop('Error: Occupancy parameter is larger than space parameter')
  if (prob < 0)                              stop('Error: Probability parameter must be between zero and one')
  if (prob > 1)                              stop('Error: Probability parameter must be between zero and one')
  
  #Compute for special case where k = 0
  if (k == 0) {
    OUT <- rep(0, n)
    return(OUT) }

  #Generate random values
  GEOM <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    if (m == Inf) { GEOM[, i] <- rgeom(n, prob = prob) }
    if (m < Inf)  { GEOM[, i] <- rgeom(n, prob = prob*(m-i+1)/m) } }
  OUT <- rowSums(GEOM)

  #Give output
  OUT }
