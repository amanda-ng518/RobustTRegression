#' Run all \eqn{\nu} estimation methods on a single simulated dataset
#'
#' @description
#' This function applies five different approaches for estimating the
#' degrees of freedom parameter \eqn{\nu} in a Student-t regression model
#' to a single simulated dataset. Each approach represents a distinct
#' estimation strategy based on profile likelihood, adjusted profile
#' likelihood, integrated Jeffreys priors, marginalized priors, or
#' joint block updates of \eqn{\nu}.
#'
#' The function provides a unified interface to compare these estimation
#' methods under the same data-generating conditions.
#'
#' @param y Numeric response vector.
#' @param x Numeric design matrix (including an intercept column).
#' @param omega_init Starting value for \code{omega} optimization.
#'   \eqn{\omega = 1 / \nu}.If simulated t-error data is used, use true 1/\eqn{\nu} for quicker convergence.
#'
#' @return A named list containing the results from all five estimation approaches:
#' \describe{
#'   \item{profile}{Result from the profile likelihood estimation.}
#'   \item{adj_profile}{Result from the adjusted profile likelihood approach.}
#'   \item{IJ}{Result from the independence Jeffrey's prior method.}
#'   \item{mar_IJ}{Result from the marginalized independence Jeffrey's prior method.}
#'   \item{nu_block}{Result from the marginalized Fisher approach for \eqn{\nu}.}
#' }
#'
#' @seealso
#' \code{\link{estimate_nu_profile}}, \code{\link{estimate_nu_adj_profile}},
#' \code{\link{estimate_nu_IJ}}, \code{\link{estimate_nu_mar_IJ}},
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
