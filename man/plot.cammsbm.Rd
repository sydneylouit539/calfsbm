% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_functions.R
\name{plot.cammsbm}
\alias{plot.cammsbm}
\title{Plot method for cammsbm objects}
\usage{
\method{plot}{cammsbm}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"cammsbm"} produced by 
\code{\link{cammsbm_vem}}.}

\item{...}{Additional arguments (currently none).}
}
\value{
Base R plot of the ELBO convergence
}
\description{
Provides a plot of the ELBO across each iteration of the VEM Algorithm
}
