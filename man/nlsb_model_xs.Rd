% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_model.R
\name{nlsb_model_xs}
\alias{nlsb_model_xs}
\title{Calculate NLSB model cross section}
\usage{
nlsb_model_xs(x, z0, a1, a2, s1, s2, wsb, x0, xf)
}
\arguments{
\item{x}{vector of x coordinates}

\item{z0}{minimum bed elevation}

\item{a1}{multiplicative parameter of the first nonlinear model}

\item{a2}{multiplicative parameter of the second nonlinear model}

\item{s1}{exponent of the first nonlinear model}

\item{s2}{exponent of the second nonlinear model}

\item{wsb}{width at the breakpoint}

\item{x0}{minimum x value, default is 0}

\item{xf}{maximum x value}
}
\value{
b modeled bed elevation
}
\description{
Calculate b(x) for the nonlinear slope break model
}
\details{
Runs in a loop, could be vectorized for speed
}
