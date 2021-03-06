% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_true_A_WP.R
\name{calc_A}
\alias{calc_A}
\title{Calculate Flow Area}
\usage{
calc_A(x, b, WSE, N = length(x), method = "trap")
}
\arguments{
\item{x}{distance from horizontal datum}

\item{b}{channel bed elevation}

\item{WSE}{water surface elevation}

\item{N}{spatial discretization, do not modify, except possibly for Riemann method}

\item{method}{"trap" "riemann" or "simpson" Default is trapezoidal.}
}
\value{
A true cross sectional flow area
}
\description{
Calculates cross sectional area of flow for a particular water surface elevation, ranging from empty channel to bankfull conditions.
}
\details{
Takes cross section info as input. 
calc_A_from_WSEw() is more efficient if WSE-w info is available
Uses a numerical method to approximate the integral of the difference between WSE and bed elevation
The numerical method only counts area where the WSE is higher than the bed elevation
Needs more work to figure out how to implement arbitrary spatial discretization
Need to do code testing, too. Confirmed that riemann and trap methods give similar results, but are they correct?
Combine with calc_WP function/put in the same R file
}
\examples{
calc_A(x = cross_sections$x[[1]], b = cross_sections$b[[1]], WSE = 138)
x <- c(1,1,1,2,3,4,4,4)
b <- c(3,2,1,1,1,1,2,3)
A <- calc_A(x, b, WSE = 3)
A <- calc_A(x, b, WSE = 2)
}
