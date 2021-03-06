% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_prediction_error.R
\name{calc_depth_prediction_error}
\alias{calc_depth_prediction_error}
\title{Calculate depth prediction error}
\usage{
calc_depth_prediction_error(pred, truth, h = NULL, type)
}
\arguments{
\item{pred}{minimum bed elevation z0 predictions}

\item{truth}{true minimum bed elevation}

\item{h}{bankfull water surface elevation. For UMESC data, it is the reference WSE. Not needed for absolute error.}

\item{type}{can be "absolute" or "relative" or "mersel"}
}
\description{
Calculates depth prediction error, either absolute or as a fraction of bankfull depth
}
\details{
Option "mersel" calculates standard error relative to bankfull depth, see Mersel et al. (2013)
}
\examples{
refWSE <- 470*0.3048 # meters, for pool 21
z0.l.error.abs <- calc_depth_prediction_error(z0.l, z0.true.ra, type = "absolute")
z0.l.error.rel <- calc_depth_prediction_error(z0.l, z0.true.ra, h = refWSE, type = "relative")
load("./Outputs/Final/p21/z0_pred.rda")
z0.true.ra <- readRDS("./Outputs/Final/p21/z0_true_ra.rds")
truth <- z0.true.ra
pred <- z0.l
z0.L.SE <- calc_depth_prediction_error(z0.l, z0.true.ra, h = refWSE, type = "mersel")
}
