#' Generates simulations from the extended balls-in-bins process
#'
#' \code{rballbin} generates simulated data from the extended balls-in-bins process.
#'
#' This function generates a simulated set of data from the extended balls-in-bins process.  The outcome is an
#' object of class \code{ballbin} with its own customised print and plot methods.  The object contains outcomes
#' of each simulation, including its effective sample-size, occupancy number, max-count number, and hitting times
#' for each possible occupancy number.
#'
#' @usage \code{rballbin(n, size, space, prob)}
#' @param n The number of simulations of the process
#' @param size The size parameter for the occupancy distribution (number of balls)
#' @param space The space pararmeter for the occupancy distribution (number of bins)
#' @param prob The probability parameter for the occupancy distribution (probability of ball occupying its bin)
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be an object of class \code{ballbin} containing \code{n} random samples from the process.

rballbin <- function(n, size, space, prob) {

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
  obs <- n
  n <- as.integer(size)
  m <- as.integer(space)

  #Check that parameters are in allowable range
  if (size != n)                            stop('Error: Size parameter is not an integer')
  if (n <= 0)                               stop('Error: Size parameter is negative')
  if (space != m)                           stop('Error: Space parameter is not an integer')
  if (m <= 0)                               stop('Error: Space parameter is negative')
  if ((prob < 0)|(prob > 1))                stop('Error: Probability parameter should be between zero and one')

  #Check the number of observations
  if (as.integer(n) != n)                   stop('Error: Proposed sample size is not an integer')
  if (as.integer(n) < 1)                    stop('Error: Proposed sample size is less than one')

  #Generate allocations
  ALLOCATION           <- matrix(0, nrow = obs, ncol = n)
  rownames(ALLOCATION) <- sprintf('Sample[%s]', 1:obs)
  colnames(ALLOCATION) <- sprintf('Ball[%s]', 1:n)
  for (i in 1:obs) {
    ALLOCATION[i,] <- sample.int(n = m, size = n, replace = TRUE) }

  #Generate samples
  SAMPLE           <- matrix(0, nrow = obs, ncol = n)
  rownames(SAMPLE) <- sprintf('Sample[%s]', 1:obs)
  colnames(SAMPLE) <- sprintf('Ball[%s]', 1:n)
  for (i in 1:obs) {
    SAMPLE[i,] <- ALLOCATION[i,]
    OCC        <- (runif(n) < prob)
    SAMPLE[i, !OCC] <- NA }

  #Generate counts and effective sample size
  COUNTS           <- matrix(0, nrow = obs, ncol = m)
  rownames(COUNTS) <- sprintf('Sample[%s]', 1:obs)
  colnames(COUNTS) <- sprintf('Bin[%s]', 1:m)
  for (i in 1:obs) {
  for (k in 1:m) {
    EFF <- (!is.na(SAMPLE[i, ]))
    COUNTS[i, k]   <- sum(SAMPLE[i, EFF] == k) } }

  #Generate effective size, occupancy and max-count
  EFF.SIZE         <- rep(0, obs)
  OCC              <- rep(0, obs)
  MAXCOUNT         <- rep(0, obs)
  names(EFF.SIZE)  <- sprintf('Sample[%s]', 1:obs)
  names(OCC)       <- sprintf('Sample[%s]', 1:obs)
  names(MAXCOUNT)  <- sprintf('Sample[%s]', 1:obs)
  for (i in 1:obs) {
    EFF.SIZE[i] <- sum(COUNTS[i, ])
    OCC[i]      <- sum(COUNTS[i, ] > 0)
    MAXCOUNT[i]      <- max(COUNTS[i, ]) }

  #Generate excess hitting times
  HITTING           <- matrix(NA, nrow = obs, ncol = m)
  rownames(HITTING) <- sprintf('Sample[%s]', 1:obs)
  colnames(HITTING) <- sprintf('Occ[%s]', 1:m)
  for (i in 1:obs) {
    tt  <- 0
    occ <- 0
    OCCINDS <- rep(FALSE, m)
    while ((occ < min(n,m))&(tt < n)) {
      tt  <- tt+1
      VAL <- SAMPLE[i, tt]
      if (!is.na(VAL)) {
      if (!OCCINDS[VAL]) {
        OCCINDS[VAL] <- TRUE
        occ <- occ+1
        HITTING[i, occ] <- tt } } } }
  for (k in 1:m) {
    HITTING[, k] <- HITTING[, k] - k }

  #Create output list
  OUT <- list(Allocation = ALLOCATION, Sample = SAMPLE, Counts = COUNTS,
              Eff.Size = EFF.SIZE, Occupancy = OCC, MaxCount = MAXCOUNT,
              Excess.Hitting = HITTING)
  class(OUT) <- 'ballbin'
  attr(OUT, 'size')  <- n
  attr(OUT, 'space') <- m
  attr(OUT, 'prob')  <- prob

  #Return output
  OUT }


print.ballbin <- function(object) {

  #Check input class
  if (!('ballbin' %in% class(object)))      stop('Error: This print method is for \'ballbin\' objects')

  #Extract information
  n         <- attributes(object)$size
  m         <- attributes(object)$space
  prob      <- attributes(object)$prob
  COUNTS    <- object$Counts
  EFF.SIZE  <- as.matrix(object$Eff.Size,  ncol = 1, drop = FALSE)
  OCCUPANCY <- as.matrix(object$Occupancy, ncol = 1, drop = FALSE)
  MAXCOUNT  <- as.matrix(object$MaxCount,  ncol = 1, drop = FALSE)
  HITTING   <- object$Excess.Hitting
  colnames(EFF.SIZE)  <- 'Eff.Size'
  colnames(OCCUPANCY) <- 'Occupancy'
  colnames(MAXCOUNT)  <- 'MaxCount'
  obs       <- nrow(COUNTS)

  #Set print parameters
  L1 <- paste0(rep('-', floor(3.5*min(max(0, m-3), 12) + 0.5*min(max(0, m-9), 6))),   collapse = '')
  L2 <- paste0(rep('-', ceiling(3.5*min(max(0, m-3), 12) + 0.5*min(max(0, m-9), 6))), collapse = '')

  #Print heading
  cat('\n      Extended Balls-in-Bins Process \n \n')
  cat(paste0('Process with ', n, ' balls randomly allocated to ', m,
             ' bins, with occupancy probability = ', round(prob, 4), ' \n \n \n'))

  #Print output for counts, etc.
  cat(paste0(L1, '-------Counts, effective sample size, occupancy and max-count------', L2), '\n \n')
  DASHES <- matrix(rep('|', obs), nrow = obs, dimnames = list(NULL, '|'))
  print(cbind(DASHES, COUNTS, DASHES, EFF.SIZE, DASHES, OCCUPANCY, DASHES, MAXCOUNT), quote = FALSE)
  cat('\n')

  #Print output for excess hitting times
  HIT.PRINT <- matrix(as.character(HITTING), nrow = obs)
  rownames(HIT.PRINT) <- rownames(HITTING)
  colnames(HIT.PRINT) <- colnames(HITTING)
  for (i in 1:obs) {
    ALL.HITS <- na.omit(HITTING[i,])
    HH <- length(ALL.HITS)
    if (HH == 0) {
      HIT.PRINT[i, ] <- rep(paste0(n, '+'), m) }
    if ((HH > 0)&(HH < m))  {
      HIT.PRINT[i, (HH+1):m] <- rep(paste0(n-HH, '+'), m-HH) } }
  cat(paste0(L1, '------------------------Excess hitting times-----------------------', L2), '\n \n')
  print(HIT.PRINT, quote = FALSE, right = TRUE)
  cat('\n')
  if (sum(is.na(HITTING)) > 0) {
    cat('(Values ending in + are censored values ---i.e., excess hitting times at least as large as the stated values.) \n \n') } }


plot.ballbin <- function(object) {

  #Check inputs
  if (!('ballbin' %in% class(object)))                                 { stop('Error: This plot method is for ballbin objects') }

  #Extract information
  n         <- attributes(object)$size
  m         <- attributes(object)$space
  prob      <- attributes(object)$prob
  COUNTS    <- object$Counts
  EFF.SIZE  <- as.matrix(object$Eff.Size,  ncol = 1, drop = FALSE)
  OCCUPANCY <- as.matrix(object$Occupancy, ncol = 1, drop = FALSE)
  MAXCOUNT  <- as.matrix(object$MaxCount,  ncol = 1, drop = FALSE)
  HITTING   <- object$Excess.Hitting
  colnames(EFF.SIZE)  <- 'Eff.Size'
  colnames(OCCUPANCY) <- 'Occupancy'
  colnames(MAXCOUNT)  <- 'MaxCount'
  obs       <- nrow(COUNTS)

  #Check installed packages and load them
  GGPLOT2 <- requireNamespace('ggplot2', quietly = TRUE)
  GREXTRA <- requireNamespace('gridExtra', quietly = TRUE)
  if (GGPLOT2) { library(ggplot2)   } else { stop('Error: Plotting a ballbins object requires the ggplot2 package') }
  if (GREXTRA) { library(gridExtra) } else { stop('Error: Plotting a ballbins object requires the gridExtra package') }

  #Create the empirical occupancy plot
  PLOTDATA1 <- data.frame(Occupancy = 0:min(n,m), Probability = 0)
  for (k in 0:min(n,m)) { PLOTDATA1$Probability[k+1] <- sum(OCCUPANCY == k)/obs }
  PLOT1 <- ggplot2::ggplot(ggplot2::aes(x = Occupancy, y = Probability, fill = 'Blue'), data = PLOTDATA1) +
           ggplot2::geom_bar() +
           ggplot2::ggtitle('Empirical Occupancy Distribution')

  #Create the max-count plot
  PLOTDATA2 <- data.frame(MaxCount = 0:n, Probability = 0)
  for (nn in 0:n) { PLOTDATA2$Probability[nn+1] <- sum(MAXCOUNT == nn)/obs }
  PLOT2 <- ggplot2::ggplot(ggplot2::aes(x = MaxCount, y = Probability, fill = 'Red'), data = PLOTDATA2) +
           ggplot2::geom_bar() +
           ggplot2::ggtitle('Empirical Max-Count Distribution') +
           ggplot2::xlab('Maximum Count')

  #Create joint plot
  SUBTITLE <- paste0('(Using ', obs, ' simulations of process with ', n, ' balls in ', m,
                     ' bins, with occupancy probability = ', round(prob, 4), ')')
  FULLPLOT <- gridExtra::grid.arrange(PLOT1, PLOT2, top = 'Balls-in-Bins Plots', bottom = SUBTITLE)

  #Print the plot
  plot(FULLPLOT) }
