% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_functions.R
\name{find_elbo}
\alias{find_elbo}
\title{Function to find the ELBO for the CAMM-SBM model.}
\usage{
find_elbo(beta_vec, network, gamma, cluster_intercept = TRUE)
}
\arguments{
\item{beta_vec}{beta_0 an beta parameters flattened into 
a combined vector}

\item{network}{List object containing adjacency (A), similarity (dis) and 
covariates (X)}

\item{gamma}{Matrix of clustering probabilities found in the E-Step}

\item{cluster_intercept}{Whether to have the intercept
term represented as a matrix (default = \code{TRUE})}
}
\value{
Observed ELBO with the given parameters
}
\description{
Function to find the ELBO for the CAMM-SBM model.
}
