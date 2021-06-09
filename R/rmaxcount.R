#' @rdname dmaxcount.all
rmaxcount <- function(n, size, space, prob = 1) {

  #Check that argument and parameters are appropriate type
  if (!is.numeric(n))                       stop('Error: Argument n is not numeric')
  if (!is.numeric(size))                    stop('Error: Size parameter is not numeric')
  if (!is.numeric(space))                   stop('Error: Space parameter is not numeric')
  if (!is.numeric(prob))                    stop('Error: Probability parameter is not numeric')

  #Check that parameters are atomic
  if (length(size)  != 1)                   stop('Error: Size parameter should be a single number')
  if (length(space) != 1)                   stop('Error: Space parameter should be a single number')
  if (length(prob)  != 1)                   stop('Error: Probability parameter should be a single number')

  #Set parameters
  nn <- as.integer(size)
  if (space == Inf) { m <- Inf } else { m <- as.integer(space) }

  #Check that parameters are in allowable range
  if (length(n) != 1)                       stop('Error: Argument n should be a single positive integer')
  if (as.integer(n) != n)                   stop('Error: Argument n should be a positive integer')
  if (min(n) < 1)                           stop('Error: Argument n should be a positive integer')
  if (size != nn)                           stop('Error: Size parameter is not an integer')
  if (nn < 0)                               stop('Error: Size parameter should be non-negative')
  if (space != m)                           stop('Error: Space parameter is not an integer')
  if (m <= 0)                               stop('Error: Space parameter should be positive')
  if ((prob < 0)|(prob > 1))                stop('Error: Probability parameter is not between zero and one')

  #Deal with trivial case where size = 0
  if (nn == 0) {
    OUT <- rep(0, n)
    return(OUT) }

  #Generate random values
  OUT <- numeric(n)
  for (i in 1:n) {
    #Generate outcomes of random sampling from extended occupancy model
    BALLS <- sample.int(m, size = nn, replace = TRUE)
    STAY  <- (runif(nn) <= prob)
    BALLS[!STAY] <- NA
    COUNTS <- rep(0, m)
    for (k in 1:m) { COUNTS[k] <- sum(BALLS == k) }
    OUT[i] <- max(COUNTS) }

  #Return the output
  OUT }
