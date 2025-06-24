#' Print method for cammsbm objects
#'
#' Displays a concise summary of a fitted CAMM-SBM model.
#'
#' @param x An object of class \code{"cammsbm"} returned by 
#' \code{\link{cammsbm_vem}}.
#' @param ... Additional arguments (currently none).
#'
#' @return Invisibly returns the input object.
#' @export
#' @method print cammsbm
print.cammsbm <- function(x, ...) {
  cat("CAMM-SBM VEM Fit\n")
  cat("---------------\n")
  cat("Final ELBO:   ", round(x$elbo, 4), "\n")
  cat("Iterations:   ", x$iter, "\n")
  cat("Time (sec):   ", round(x$time, 2), "\n")
  invisible(x)
}


#' Summary method for cammsbm objects
#'
#' Provides a detailed summary of the estimated parameters and model fit.
#'
#' @param object An object of class \code{"cammsbm"} produced by \code{\link{cammsbm_vem}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#' @export
#' @method summary cammsbm
summary.cammsbm <- function(object, ...) {
  cat("Summary of CAMM-SBM VEM fit\n")
  cat("---------------------------\n")
  cat("Number of clusters (K): ", ncol(object$gamma), "\n")
  cat("Iterations until convergence: ", object$iter, "\n")
  cat("Final ELBO: ", round(object$elbo, 4), "\n")
  cat("Run time (sec): ", round(object$time, 2), "\n\n")
  
  cat("Beta0 Estimates:\n")
  print(object$beta0)
  cat("\nBeta Estimates:\n")
  print(object$beta)
  
  invisible(object)
}


#' Plot method for cammsbm objects
#'
#' Provides a plot of the ELBO across each iteration of the VEM Algorithm
#'
#' @param x An object of class \code{"cammsbm"} produced by 
#' \code{\link{cammsbm_vem}}.
#' @param ... Additional arguments (currently none).
#'
#' @return Base R plot of the ELBO convergence
#' @export
#' @method plot cammsbm
plot.cammsbm <- function(x, ...) {
  plot(x$evals, pch = 19,
       xlab = "Iteration", ylab = "ELBO",
       main = "Plot of ELBO Progression")
}






