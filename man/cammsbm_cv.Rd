% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_functions.R
\name{cammsbm_cv}
\alias{cammsbm_cv}
\title{Function to perform cross-validation on the CAMM-SBM model. This function can  
be used to select the number of clusters or to perform feature selection}
\usage{
cammsbm_cv(network, K, folds = 5, mcmc = FALSE)
}
\arguments{
\item{network}{List object containing adjacency (A), similarity (dis), 
and covariates (X)}

\item{K}{List object containing adjacency (A), similarity (dis), 
and covariates (X)}

\item{folds}{Number of folds of cross-validation (default = 5)}

\item{mcmc}{Boolean indicating whether to run the model in MCMC format.
The model will run much slower, but may be more accurate}
}
\value{
Combined BIC of all the test sets
}
\description{
Function to perform cross-validation on the CAMM-SBM model. This function can  
be used to select the number of clusters or to perform feature selection
}
