#' Generates simulations from the extended balls-in-bins process
#'
#' \code{sample.ballbin} generates simulated data from the extended balls-in-bins process.
#'
#' This function generates a simulated set of data from the extended balls-in-bins process.  The outcome is an
#' object of class \code{ballbin} containing the simulations from the process.  The output object contains the
#' initial bin-allocation and resulting samples from the process.  Calling \code{summary} on the simulation object
#' creates a new object of class \code{summary.ballbin} containing summary statistics for each sample, including
#' the bin-counts, effective sample-size, occupancy number, max-count number, and hitting times for each possible
#' occupancy number.  Each of these objects has a custom printing and plot methods to give user-friendly output.
#'
#' @param n The number of simulations of the process
#' @param size The size parameter for the occupancy distribution (number of balls)
#' @param space The space pararmeter for the occupancy distribution (number of bins)
#' @param prob The probability parameter for the occupancy distribution (probability of ball occupying its bin)
#' @param alloc.prob (Optional) A probability vector for the allocation probabilities for the bins
#' @return If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
#' be a list of class \code{ballbin} containing \code{n} random samples from the process.  If you call \code{summary}
#' on this object the output will be another list of class \code{summary.ballbin} containing summary statistics
#' for each random sample from the process.

sample.ballbin <- function(n, size, space, prob, alloc.prob = NULL) {

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

  #Check the allocation probability vector
  if (!is.null(alloc.prob)) {
    if (!is.numeric(alloc.prob))            stop('Error: Allocation probability vector is not numeric')
    if (length(alloc.prob) != m)            stop('Error: Allocation probability vector should have length equal to the space parameter')
    if (min(alloc.prob) < 0)                stop('Error: Allocation probability vector has one or more negative elements')
    if (sum(alloc.prob) != 1)               stop('Error: Allocation probability vector does not sum to one') }

  #Generate samples
  ALLOCATION       <- matrix(0, nrow = obs, ncol = n)
  rownames(ALLOCATION) <- sprintf('Sample[%s]', 1:obs)
  colnames(ALLOCATION) <- sprintf('Ball[%s]', 1:n)
  for (i in 1:obs) {
    if (is.null(alloc.prob)) {
      ALLOCATION[i,] <- sample.int(n = m, size = n, replace = TRUE)
    } else {
      ALLOCATION[i,] <- sample.int(n = m, size = n, replace = TRUE, prob = alloc.prob) } }

  #Generate samples
  SAMPLE           <- matrix(0, nrow = obs, ncol = n)
  rownames(SAMPLE) <- sprintf('Sample[%s]', 1:obs)
  colnames(SAMPLE) <- sprintf('Ball[%s]', 1:n)
  for (i in 1:obs) {
    SAMPLE[i,] <- ALLOCATION[i,]
    OCC        <- (runif(n) < prob)
    SAMPLE[i, !OCC] <- NA }

  #Create output list
  if (is.null(alloc.prob)) {
    PARS <- list(size = n, space = m, prob = prob)
  } else {
    PARS <- list(size = n, space = m, prob = prob, alloc.prob = alloc.prob) }
  OUT  <- list(parameters = PARS, allocation = ALLOCATION, sample = SAMPLE)
  class(OUT) <- 'ballbin'

  #Return output
  OUT }

#' @describeIn sample.ballbin prints the sample
print.ballbin <- function(x, ...) {

  #Check input class
  if (!('ballbin' %in% class(x)))      stop('Error: This print method is for \'ballbin\' objects')

  #Extract information
  PARS       <- x$parameters
  n          <- PARS$size
  m          <- PARS$space
  prob       <- PARS$prob
  alloc.prob <- PARS$alloc.prob
  SAMPLE     <- x$sample
  obs        <- nrow(SAMPLE)

  #Print heading
  if (prob == 1) {
    cat('\n      Balls-in-Bins Process \n \n')
    cat(paste0('Process with ', n, ' balls randomly allocated to ', m, ' bins \n')) }
  if (prob < 1) {
    cat('\n      Extended Balls-in-Bins Process \n \n')
    cat(paste0('Process with ', n, ' balls randomly allocated to ', m,
               ' bins, with occupancy probability = ', round(prob, 4), ' \n')) }
  if (!is.null(alloc.prob)) {
    PP <- round(alloc.prob, 4)
    cat('Allocation probability vector = (')
    cat(PP, sep = ', ')
    cat(') \n') }
  cat('\n \n')

  #Print output
  print(SAMPLE)
  cat('\n')}

#' @describeIn sample.ballbin plots the sample
plot.ballbin <- function(x, ..., ball.size = NULL, ball.color = NULL, ball.colour = ball.color, max.plots = 30) {

  #Check inputs
  if (!('ballbin' %in% class(x)))            stop('Error: This plot method is for \'ballbin\' objects')
  if (!is.null(ball.size)) {
    if (!is.numeric(ball.size))                   stop('Error: ball.size must be a positive number')
    if (length(ball.size) != 1)                   stop('Error: ball.size must be a single positive number')
    if (min(ball.size) <= 0)                      stop('Error: ball.size must be a positive number') }
  if ((!missing(ball.color))&(!missing(ball.colour))) {
    stop('Error: Specify ball.color or ball.colour but not both') }
  if (!is.null(ball.colour)) {
    if (!is.character(ball.colour))               stop('Error: ball.color must be in \'colours()\'')
    if (length(ball.colour) != 1)                 stop('Error: ball.color must be in \'colours()\'')
    if (!(ball.colour %in% colours()))            stop('Error: ball.color must be in \'colours()\'') }
  if (!is.numeric(max.plots))                     stop('Error: max.plots must be a positive integer')
  if (length(max.plots) != 1)                     stop('Error: max.plots must be a single positive integer')
  MAX.PLOTS <- as.integer(max.plots)
  if (max.plots != MAX.PLOTS)                     stop('Error: max.plots must be a positive integer')

  #Extract information
  PARS       <- x$parameters
  n          <- PARS$size
  m          <- PARS$space
  prob       <- PARS$prob
  alloc.prob <- PARS$alloc.prob
  ALLOC      <- x$allocation
  SAMPLE     <- x$sample
  obs        <- nrow(SAMPLE)

  #Limit plot size to 100 plots
  if (obs > MAX.PLOTS) {
    warning(paste0('There were more than ', MAX.PLOTS, ' samples --- plot shows only the first ', MAX.PLOTS, ' samples'))
    ALLOC  <- ALLOC[1:MAX.PLOTS, ]
    SAMPLE <- SAMPLE[1:MAX.PLOTS, ]
    obs    <- MAX.PLOTS }

  #Check installed packages and load them
  GGPLOT2 <- requireNamespace('ggplot2', quietly = TRUE)
  GREXTRA <- requireNamespace('gridExtra', quietly = TRUE)
  if (GGPLOT2) { library(ggplot2)   } else { stop('Error: Plotting a \'ballbins\' object requires the ggplot2 package') }
  if (GREXTRA) { library(gridExtra) } else { stop('Error: Plotting a \'ballbins\' object requires the gridExtra package') }

  #Generate occupancy indicators
  OCC           <- matrix(FALSE, nrow = obs, ncol = m)
  rownames(OCC) <- sprintf('Sample[%s]', 1:obs)
  colnames(OCC) <- sprintf('Bin[%s]', 1:m)
  for (i in 1:obs) {
    for (k in 1:m) {
      EFF <- (!is.na(SAMPLE[i, ]))
      OCC[i, k]   <- (sum(SAMPLE[i, EFF] == k) > 0) } }

  #Set ball and bin size and colour
  if (!is.null(ball.size))   { BALL.SIZE   <- ball.size   } else { BALL.SIZE   <- 2 }
  if (!is.null(ball.colour)) { BALL.COLOUR <- ball.colour } else { BALL.COLOUR <- 'blue' }

  #Set subtitle for plot
  if (obs == 1) {
    SUBTITLE <- paste0('(Showing one simulation of process with ', n, ' balls allocated to ', m,
                       ' bins, with occupancy probability = ', round(prob, 4), ')') }
  if (obs > 1) {
    SUBTITLE <- paste0('(Showing ', obs, ' simulations of process with ', n, ' balls allocated to ', m,
                       ' bins, with occupancy probability = ', round(prob, 4), ')') }
  if (!is.null(alloc.prob)) {
    SUBTITLE <- paste0(SUBTITLE, '\n', '(Allocation probability vector = ', round(alloc.prob, 4), ')') }
  if (prob < 1) {
    SUBTITLE <- paste0(SUBTITLE, '\n', '(Filled balls are those that occupy their allocated bins --- filled bins are those that are occupied)') }

  #Create the ball data
  BALLDATA <- data.frame(Sample = rep(1:obs, each = n), Ball = rep(1:n, times = obs), Bin = 0, Occ = TRUE, Type = 'Ball',
                         Label  = rep(sprintf('Sample[%s]', 1:obs), each = n))
  for (i in 1:(obs*n)) {
    BALLDATA$Bin[i] <- ALLOC[BALLDATA$Sample[i], BALLDATA$Ball[i]]
    BALLDATA$Occ[i] <- !is.na(SAMPLE[BALLDATA$Sample[i], BALLDATA$Ball[i]]) }
  BALLDATA$Shape <- factor(rep('21', obs*n), levels = c('21', '22'))

  #Create the bin data
  BINDATA <- data.frame(Sample = rep(1:obs, each = m), Ball = n+1+(BALL.SIZE)/4, Bin = rep(1:m, times = obs), Occ = c(t(OCC)), Type = 'Bin',
                        Label  = rep(sprintf('Sample[%s]', 1:obs), each = m))
  BINDATA$Shape  <- factor(rep('22', obs*m), levels = c('21', '22'))

  #Create the sample plots
  PLOTDATA <- rbind(BALLDATA, BINDATA)
  LABELS <- as.list(rownames(ALLOC))
  names(LABELS) <- 1:obs
  PLOT <- ggplot2::ggplot(ggplot2::aes(x = Bin, y = Ball, fill = Occ), data = PLOTDATA) +
          ggplot2::geom_point(size = BALL.SIZE,   shape = 21, data = BALLDATA) +
          ggplot2::geom_point(size = BALL.SIZE+2, shape = 22, data = BINDATA) +
          ggplot2::scale_x_continuous(labels = 1:m, breaks = 1:m) +
          ggplot2::scale_y_reverse() +
          ggplot2::facet_wrap( ~ factor(Label, levels = LABELS)) +
          ggplot2::scale_fill_manual(values = c('#FFFFFF00', BALL.COLOUR), guide = 'none') +
          ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold'),
                         panel.grid.minor.x = ggplot2::element_blank(),
                         axis.title.x  = ggplot2::element_blank(),
                         axis.text.x   = ggplot2::element_blank(),
                         axis.ticks.x  = ggplot2::element_blank(),
                         axis.title.y  = ggplot2::element_blank(),
                         axis.text.y   = ggplot2::element_blank(),
                         axis.ticks.y  = ggplot2::element_blank(),
                         panel.grid.major.y = ggplot2::element_blank(),
                         panel.grid.minor.y = ggplot2::element_blank(),
                         panel.spacing = unit(1, 'lines')) +
          ggplot2::ggtitle('Balls-in-Bins Sample Plots') +
          ggplot2::labs(subtitle = SUBTITLE)

  #Print the plot
  plot(PLOT) }

#' @describeIn sample.ballbin summarizes the sample
summary.ballbin <- function(object, ...) {

  #Check input class
  if (!('ballbin' %in% class(object)))      stop('Error: This summary method is for \'ballbin\' objects')

  #Extract information
  PARS   <- object$parameters
  n      <- PARS$size
  m      <- PARS$space
  SAMPLE <- object$sample
  obs    <- nrow(SAMPLE)

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
  OUT <- list(parameters = PARS, sample = SAMPLE, counts = COUNTS, eff.size = EFF.SIZE,
              occupancy = OCC, max.count = MAXCOUNT, excess.hitting = HITTING)
  class(OUT) <- 'summary.ballbin'

  #Return output
  OUT }

#' @describeIn sample.ballbin prints the summary
print.summary.ballbin <- function(x, ...) {

  #Check input class
  if (!('summary.ballbin' %in% class(x)))      stop('Error: This print method is for \'summary.ballbin\' objects')

  #Extract information
  PARS       <- x$parameters
  n          <- PARS$size
  m          <- PARS$space
  prob       <- PARS$prob
  alloc.prob <- PARS$alloc.prob
  COUNTS     <- x$counts
  EFF.SIZE   <- as.matrix(x$eff.size,  ncol = 1, drop = FALSE)
  OCCUPANCY  <- as.matrix(x$occupancy, ncol = 1, drop = FALSE)
  MAXCOUNT   <- as.matrix(x$max.count, ncol = 1, drop = FALSE)
  HITTING    <- x$excess.hitting
  colnames(EFF.SIZE)  <- 'Eff.Size'
  colnames(OCCUPANCY) <- 'Occupancy'
  colnames(MAXCOUNT)  <- 'MaxCount'
  obs       <- nrow(COUNTS)

  #Set print parameters
  L1 <- paste0(rep('-', floor(3.5*min(max(0, m-3), 12) + 0.5*min(max(0, m-9), 6))),   collapse = '')
  L2 <- paste0(rep('-', ceiling(3.5*min(max(0, m-3), 12) + 0.5*min(max(0, m-9), 6))), collapse = '')

  #Print heading
  if (prob == 1) {
    cat('\n      Balls-in-Bins Process \n \n')
    cat(paste0('Process with ', n, ' balls randomly allocated to ', m, ' bins \n')) }
  if (prob < 1) {
    cat('\n      Extended Balls-in-Bins Process \n \n')
    cat(paste0('Process with ', n, ' balls randomly allocated to ', m,
               ' bins, with occupancy probability = ', round(prob, 4), ' \n')) }
  if (!is.null(alloc.prob)) {
    PP <- round(alloc.prob, 4)
    cat('Allocation probability vector = (')
    cat(PP, sep = ', ')
    cat(') \n') }
  cat('\n \n')

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
    cat('(Values ending in + are right-censored values.) \n \n') } }

#' @describeIn sample.ballbin plots the summary
#' @param x,bar.color,bar.colour plotting arguments
plot.summary.ballbin <- function(x, ..., bar.color = NULL, bar.colour = bar.color) {

  #Check inputs
  if (!('summary.ballbin' %in% class(x)))              stop('Error: This plot method is for \'summary.ballbin\' objects')
  if ((!missing(bar.color))&(!missing(bar.colour))) {
    stop('Error: Specify bar.color or bar.colour but not both') }
  if (!is.null(bar.colour)) {
    if (!is.character(bar.colour))                          stop('Error: bar.color must be in \'colours()\'')
    if (length(bar.colour) != 1)                            stop('Error: bar.color must be in \'colours()\'')
    if (!(bar.colour %in% colours()))                       stop('Error: bar.color must be in \'colours()\'') }

  #Extract information
  PARS       <- x$parameters
  n          <- PARS$size
  m          <- PARS$space
  prob       <- PARS$prob
  alloc.prob <- PARS$alloc.prob
  COUNTS     <- x$counts
  EFF.SIZE   <- as.matrix(x$eff.size,  ncol = 1, drop = FALSE)
  OCCUPANCY  <- as.matrix(x$occupancy, ncol = 1, drop = FALSE)
  MAXCOUNT   <- as.matrix(x$max.count, ncol = 1, drop = FALSE)
  HITTING    <- x$excess.hitting
  colnames(EFF.SIZE)  <- 'Eff.Size'
  colnames(OCCUPANCY) <- 'Occupancy'
  colnames(MAXCOUNT)  <- 'MaxCount'
  obs       <- nrow(COUNTS)

  #Check installed packages and load them
  GGPLOT2 <- requireNamespace('ggplot2', quietly = TRUE)
  GREXTRA <- requireNamespace('gridExtra', quietly = TRUE)
  if (GGPLOT2) { library(ggplot2)   } else { stop('Error: Plotting a \'summary.ballbins\' object requires the ggplot2 package') }
  if (GREXTRA) { library(gridExtra) } else { stop('Error: Plotting a \'summary.ballbins\' object requires the gridExtra package') }

  #Create plot data
  PLOTDATA <- data.frame(Occupancy = 0:min(n,m), Probability = 0)
  for (k in 0:min(n,m)) { PLOTDATA$Probability[k+1] <- sum(OCCUPANCY == k)/obs }

  #Set subtitle for plot
  SUBTITLE <- paste0('(', n, ' balls allocated to ', m, ' bins, with occupancy probability = ', round(prob, 4), ')')
  if (!is.null(alloc.prob)) {
    SUBTITLE <- paste0(SUBTITLE, '\n', '(Allocation probability vector = ', round(alloc.prob, 4), ')') }

  #Create the empirical occupancy plot
  if (!is.null(bar.colour)) { BAR.COLOUR <- bar.colour } else { BAR.COLOUR <- 'blue' }
  PLOT <- ggplot2::ggplot(ggplot2::aes(x = Occupancy, y = Probability), data = PLOTDATA) +
          ggplot2::geom_bar(stat = 'identity', fill = BAR.COLOUR) +
          ggplot2::scale_x_continuous(labels = 0:m, breaks = 0:m) +
          ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold'),
                         axis.title.x  = ggplot2::element_text(hjust = 0.5, face = 'bold'),
                         axis.title.y  = ggplot2::element_text(hjust = 0.5, face = 'bold',
                                                               margin = ggplot2::margin(t = 0, r = 8, b = 0, l = 0))) +
          ggplot2::ggtitle('Empirical Occupancy Distribution') +
          ggplot2::labs(subtitle = SUBTITLE)

  #Print the plot
  plot(PLOT) }

