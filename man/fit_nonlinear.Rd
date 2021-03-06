% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_nonlinear.R
\name{fit_nonlinear}
\alias{fit_nonlinear}
\title{Nonlinear Fit for WSE-w Relationship}
\usage{
fit_nonlinear(WSEw, h = 5, maxiter = 100)
}
\arguments{
\item{WSEw}{WSEw data (at a given level of exposure)}

\item{h}{minimum amount of data to perform fit}

\item{maxiter}{maximum number of iterations for the parameter estimation algorithm}
}
\description{
Nonlinear Fit for WSE-w Relationship
}
\details{
Uses a default initial guess, but if the fit is very bad, 
attempts several fits performed with different starting guesses 
and returns the best fit.
See attributes(fit)$ef for error flags. 
Value 0 means no error,
1 means not enough data points, 
2 means singular gradient at initial guess for nlsLM
}
\examples{
WSEw <- data.frame(WSE = NA, w = NA) # test case
model1 <- fit_nonlinear(WSEw, h = 5, maxiter = 100)
}
