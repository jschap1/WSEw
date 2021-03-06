% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred_linear_par.R
\name{pred_linear_par}
\alias{pred_linear_par}
\title{Predict parameters (parallel)}
\usage{
pred_linear_par(r, exclude = FALSE)
}
\arguments{
\item{r}{cross section number}

\item{exclude}{boolean that determines whether or not to exclude predictions from physically unrealistic model fits}
}
\value{
list of z0, A, W0, and A0 predictions
}
\description{
Parallelized functions for main_d6 hydraulic parameter predictions in parallel with foreach
for WSE-w fits.
}
\details{
Loads fitted model (linear, slope break, multiple slope break, nonlinear, or nonlinear slope break)
and uses it to predict hydraulic parameters
error.flag values: 0 = no error, 1 = negative slope, 3 = no model was fit
}
