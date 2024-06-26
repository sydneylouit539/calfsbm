% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbm_functions.R
\name{gibbs_sampler_fixed_k}
\alias{gibbs_sampler_fixed_k}
\title{(Obsolete) Gibbs sampler}
\usage{
gibbs_sampler_fixed_k(K, alpha, beta0, beta, niter, A, S_ij, directed = FALSE)
}
\arguments{
\item{K}{A positive integer indicating the true number of 
clusters (should be optimized already)}

\item{alpha}{Vector of prior probabilities of being in each group}

\item{beta0}{Prior intercept for logistic regression}

\item{beta}{Prior parameters for logistic regression, in matrix form}

\item{niter}{Number of iterations to run}

\item{A}{Adjacency matrix}

\item{S_ij}{Distance matrix}

\item{directed}{logical (default = \code{FALSE}); if \code{FALSE}, the 
network is undirected}
}
\value{
List containing the estimated beta, node membership, and history
}
\description{
Gibbs sampler with fixed number of clusters
}
\note{
Deprecated
}
