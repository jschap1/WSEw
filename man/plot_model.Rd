% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_model.R
\name{plot_model}
\alias{plot_model}
\title{Plot modeled cross section}
\usage{
plot_model(model, type, WSEw = NULL, add = FALSE, ...)
}
\arguments{
\item{type}{Either "linear", "sb", "sbm", "nl", or "nlsb"}

\item{WSEw}{required for nonlinear models}

\item{add}{TRUE to add plot to current figure, FALSE to create a new figure. Default is FALSE}

\item{...}{additional parameters for plot or lines}
}
\description{
Computes modeled bed elevation and plots the cross section
}
\details{
Eventually, it might make sense for this to return the x and b values in list similar to cross_sections from auto_transects
}
\examples{
plot(cross_sections$x[[2]], cross_sections$b[[2]], type = "l")
plot_model(lf1, type = "linear", col = "red", add = TRUE)
}
