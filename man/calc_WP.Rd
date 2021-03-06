% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_true_A_WP.R
\name{calc_WP}
\alias{calc_WP}
\title{Calculate Wetted Perimeter}
\usage{
calc_WP(x, b, WSE, method = "wide")
}
\arguments{
\item{x}{distance from horizontal datum}

\item{b}{channel bed elevation}

\item{WSE}{water surface elevation}

\item{method}{"neal" for Jeff Neal's approximation in Neal et al. (2015), "wide" for wide-channel approximation}
}
\value{
WP wetted perimeter
}
\description{
Calculates wetted perimeter for a particular water surface elevation, ranging from empty channel to bankfull conditions
}
\details{
Uses a numerical method to approximate the integral of the difference between WSE and bed elevation
The numerical method only counts area where the WSE is higher than the bed elevation
Needs more work to figure out how to implement arbitrary spatial discretization
Performs poorly for a small number of data points or WSE close to bed elevation, see the example
}
\examples{
calc_WP(x = cross_sections$x[[1]], b = cross_sections$b[[1]], WSE = 138)
x <- c(1,1,1,2,3,4,4,4)
b <- c(3,2,1,1,1,1,2,3)
WP <- calc_WP(x, b, WSE = 3)
WP <- calc_WP(x, b, WSE = 2)
}
