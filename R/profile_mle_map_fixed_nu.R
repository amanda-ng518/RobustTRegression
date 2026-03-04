# ============================================================
# Helper functions for Profile MLE and MAP (fixed nu)
# ============================================================

#' Profile MLE of \code{sigma} and \code{beta} for Fixed Degrees of Freedom
#'
#' This internal helper function computes the profile maximum likelihood estimates (MLEs)
#' of the regression coefficients (\code{beta}) and the scale parameter (\code{sigma})
#' for a fixed degrees of freedom \code{nu} under a Student-\emph{t} error model.
#'
#' It returns the estimated parameters along with the corresponding negative log-likelihood
#' value evaluated at those estimates.
#'
#' @param nu Numeric. The fixed degrees of freedom for the t-distribution.
#' @param y Numeric vector. The response variable.
#' @param x Numeric matrix. The design matrix (must include an intercept column).
#'
#' @return A list containing:
#' \describe{
#'   \item{beta_hat}{Profile MLE of \code{beta}.}
#'   \item{sigma_hat}{Profile MLE of \code{sigma}.}
#'   \item{loglik}{Negative log-likelihood value evaluated at the fitted parameters.}
#'   \item{convergence}{Convergence code returned by the optimizer (0 indicates successful convergence).}
#' }
#'
#' @keywords internal
fit_profile_mle_fixed_nu <- function(y, x, nu = Inf, par_init = NULL, control = list()) {
  n <- length(y)
  p <- ncol(x)
  if (is.null(par_init)) {
    lmfit <- lm(y ~ x - 1) # x includes intercept
    sigma_init <- stats::sd(residuals(lmfit))
    par_init <- c(coef(lmfit), sigma_init)
  }
  treg_neg_loglik_no_nu <- function(par) {
    beta <- par[1:p]
    sigma <- abs(par[p+1])
    z <- (y - x %*% beta)/sigma
    if (is.infinite(nu)) {
      ll <- sum(dnorm(z, log = TRUE)) - n * log(sigma)
    } else {
      ll <- sum(dt(x = z, df = nu, log = TRUE)) - n * log(sigma)
    }
    -ll
  }
  opt <- optim(par = par_init, fn = treg_neg_loglik_no_nu, method = "BFGS", control = c(list(maxit = 10000), control))
  par_hat <- opt$par
  list(beta = par_hat[1:p], sigma = abs(par_hat[p+1]), loglik = -opt$value, convergence = opt$convergence)
}

#' Profile MAP Estimates of \code{sigma} and \code{beta} for Fixed Degrees of Freedom
#'
#' This function computes the maximum a posteriori (MAP) estimates of the regression
#' coefficients (\code{beta}) and the scale parameter (\code{sigma}) for a fixed
#' value of the degrees of freedom parameter \code{nu} in a Student-\emph{t} regression model.
#'
#' The optimization is performed over the log-posterior, which includes:
#' \itemize{
#'   \item The log-likelihood under a Student-\emph{t} error model,
#'   \item An improper prior on \code{sigma} proportional to \eqn{1/\sigma} (i.e., \eqn{-\log(\sigma)}),
#'   \item A flat prior on \code{beta} (i.e., constant log-density).
#' }
#'
#' This function is used as the first-stage optimization step in
#' estimating the full joint posterior over all model parameters, including \code{nu},
#' where the full objective includes an additional prior on \code{nu}.
#'
#' @param nu Numeric. Fixed degrees of freedom for the t-distribution.
#' @param y Numeric vector. The response variable.
#' @param x Numeric matrix. The design matrix (must include an intercept column).
#'
#' @return A list containing:
#' \describe{
#'   \item{beta_hat}{MAP estimate of the regression coefficients.}
#'   \item{sigma_hat}{MAP estimate of the scale parameter.}
#'   \item{loglik}{Value of the log-posterior objective (excluding the \code{nu} prior) evaluated at the estimates.}
#'   \item{convergence}{Convergence code from \code{\link[stats]{optim}}. A value of 0 indicates successful convergence.}
#' }
#'
#' @keywords internal
fit_profile_map_fixed_nu <- function(y, x, nu, par_init = NULL, control = list(maxit = 10000)) {
  n <- length(y)
  p <- ncol(x)

  # Default initialization via OLS
  if (is.null(par_init)) {
    lmfit <- lm(y ~ x - 1)  # x includes intercept
    beta_init <- coef(lmfit)
    sigma_init <- sd(residuals(lmfit))
    par_init <- c(beta_init, sigma_init)
  }

  # Negative log-likelihood + log sigma prior for fixed nu
  treg_neg_loglik_sigmaprior_no_nu <- function(par) {
    beta <- par[1:p]
    sigma <- abs(par[p + 1])

    z <- as.vector((y - x %*% beta) / sigma)
    if (is.infinite(nu)) {
      ll <- sum(log(dnorm(z))) - (n + 1) * log(sigma)
    } else {
      ll <- sum(log(dt(z, df = nu))) - (n + 1) * log(sigma)
    }
    -ll
  }
  opt <- optim(par = par_init, fn = treg_neg_loglik_sigmaprior_no_nu, method = "BFGS", control = c(list(maxit = 10000), control))
  par_hat <- opt$par
  list(beta = par_hat[1:p], sigma = abs(par_hat[p+1]), loglik = -opt$value, convergence = opt$convergence)
}
