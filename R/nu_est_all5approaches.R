#' Run All \eqn{\nu} Estimation Methods on a Single Simulated Dataset
#'
#' @description
#' Applies five different methods for estimating the degrees of freedom parameter
#' \eqn{\nu} in a Student-\emph{t} regression model to a single simulated dataset.
#' Each method represents a distinct estimation strategy based on:
#'
#' \itemize{
#'   \item Profile likelihood
#'   \item Adjusted profile likelihood
#'   \item Independent Jeffrey's
#'   \item Marginalized independent Jeffrey's
#'   \item Marginalized Fisher
#' }
#'
#' This function provides a unified interface to compare these approaches
#' under identical data-generating conditions.
#'
#' @param y Numeric vector. The response variable.
#' @param x Numeric matrix. The design matrix (must include an intercept column).
#' @param omega_init Numeric. Initial value for \eqn{\omega = 1 / \nu}.
#'   If using simulated data, setting this to the true value of \eqn{1 / \nu}
#'   can improve convergence.
#'
#' @return A named list containing the results from all five estimation methods:
#' \describe{
#'   \item{profile}{Result from the profile likelihood estimation.}
#'   \item{adj_profile}{Result from the adjusted profile likelihood approach.}
#'   \item{IJ}{Result from the independent Jeffrey's approach.}
#'   \item{mar_IJ}{Result from the marginalized independent Jeffrey's approach.}
#'   \item{nu_block}{Result from the marginalized Fisher approach.}
#' }
#'
#' @seealso
#' \code{\link{estimate_nu_profile}},
#' \code{\link{estimate_nu_adj_profile}},
#' \code{\link{estimate_nu_IJ}},
#' \code{\link{estimate_nu_mar_IJ}},
#' \code{\link{estimate_nu_nu_block}}
#'
#' @export
run_all_estimators <- function(y, x, omega_init = 1/2) {
  res_profile <- estimate_nu_profile(y, x, omega_init)
  res_adj <- estimate_nu_adj_profile(y, x, omega_init)
  res_IJ <- estimate_nu_IJ(y, x, omega_init)
  res_mar_IJ <- estimate_nu_mar_IJ(y, x, omega_init)
  res_nu_block <- estimate_nu_nu_block(y, x, omega_init)
  list(profile = res_profile,
       adj_profile = res_adj,
       IJ = res_IJ,
       mar_IJ = res_mar_IJ,
       nu_block = res_nu_block)
}
