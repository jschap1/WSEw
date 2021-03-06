% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_model_WP.R
\name{calc_model_WP}
\alias{calc_model_WP}
\title{Calculate modeled wetted perimeter}
\usage{
calc_model_WP(model, type, w = NULL)
}
\arguments{
\item{model}{model object}

\item{w}{flow width vector, needed for nonlinear models}
}
\value{
WP wetted perimeter
}
\description{
Calculates wetted perimeter for a modeled cross section using one of several types of geometric channel models
}
\details{
Uses correct formulas to calculate WP for the modeled cross sections. 
For all the UMRB cross sections, width is much larger than depth, so WP = w, though.
}
