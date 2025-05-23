% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbm_functions.R
\name{update_beta}
\alias{update_beta}
\title{Update Beta in Gibbs Sampler}
\usage{
update_beta(K, group, directed = FALSE, offset = FALSE)
}
\arguments{
\item{K}{A positive integer indicating the true number of clusters}

\item{group}{Model matrix to be fitted on}

\item{directed}{logical; if \code{FALSE} (default), the MCMC output is from 
an undirected network}

\item{offset}{Boolean indicating whether to use offset term}
}
\value{
beta0 and beta, with beta as a matrix
}
\description{
Helper function to update beta according to adjacency and node membership
}
\note{
Function \code{update_beta} is a helper for the initialization process 
in the \code{calf_sbm_nimble} function, deriving an estimate for beta using 
the initial clustering configuration as input
}
