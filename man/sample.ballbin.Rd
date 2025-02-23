% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.ballbin.R
\name{sample.ballbin}
\alias{sample.ballbin}
\alias{print.ballbin}
\alias{plot.ballbin}
\alias{summary.ballbin}
\alias{print.summary.ballbin}
\alias{plot.summary.ballbin}
\title{Generates simulations from the extended balls-in-bins process}
\usage{
sample.ballbin(n, size, space, prob, alloc.prob = NULL)

\method{print}{ballbin}(x, ...)

\method{plot}{ballbin}(
  x,
  ...,
  ball.size = NULL,
  ball.color = NULL,
  ball.colour = ball.color,
  max.plots = 30
)

\method{summary}{ballbin}(object, ...)

\method{print}{summary.ballbin}(x, ...)

\method{plot}{summary.ballbin}(x, ..., bar.color = NULL, bar.colour = bar.color)
}
\arguments{
\item{n}{The number of simulations of the process}

\item{size}{The size parameter for the occupancy distribution (number of balls)}

\item{space}{The space pararmeter for the occupancy distribution (number of bins)}

\item{prob}{The probability parameter for the occupancy distribution (probability of ball occupying its bin)}

\item{alloc.prob}{(Optional) A probability vector for the allocation probabilities for the bins}

\item{x, object}{ballbin objects (for generics)}

\item{...}{unused}

\item{ball.size, ball.color, ball.colour, max.plots}{Set the size, color, and number of plots}

\item{bar.color, bar.colour}{plotting arguments}
}
\value{
If all inputs are correctly specified (i.e., parameters are in allowable range) then the output will
be a list of class \code{ballbin} containing \code{n} random samples from the process.  If you call \code{summary}
on this object the output will be another list of class \code{summary.ballbin} containing summary statistics
for each random sample from the process.
}
\description{
\code{sample.ballbin} generates simulated data from the extended balls-in-bins process.
}
\details{
This function generates a simulated set of data from the extended balls-in-bins process.  The outcome is an
object of class \code{ballbin} containing the simulations from the process.  The output object contains the
initial bin-allocation and resulting samples from the process.  Calling \code{summary} on the simulation object
creates a new object of class \code{summary.ballbin} containing summary statistics for each sample, including
the bin-counts, effective sample-size, occupancy number, max-count number, and hitting times for each possible
occupancy number.  Each of these objects has a custom printing and plot methods to give user-friendly output.
}
\section{Methods (by generic)}{
\itemize{
\item \code{print}: prints the sample

\item \code{plot}: plots the sample

\item \code{summary}: summarizes the sample

\item \code{print}: prints the summary

\item \code{plot}: plots the summary
}}

\examples{

d <- sample.ballbin(12, 10, 4, .4)
print(d)
plot(d)
summary(d)
plot(summary(d))
}
