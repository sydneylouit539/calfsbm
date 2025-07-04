% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_functions.R
\name{vem_e}
\alias{vem_e}
\title{Helper E-Step for VEM algorithm}
\usage{
vem_e(network, gamma, beta0, beta)
}
\arguments{
\item{network}{List object containing adjacency (A), similarity (dis), 
and covariates (X)}

\item{gamma}{n x K matrix of clustering probabilities}

\item{beta0}{K x K intercept matrix}

\item{beta}{K x K matrix}
}
\value{
Fitted n x K matrix of clustering probabilities
}
\description{
Helper E-Step for VEM algorithm
}
