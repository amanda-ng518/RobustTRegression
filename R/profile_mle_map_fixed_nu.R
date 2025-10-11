#' Profile MLE of Sigma and Beta for Fixed Degrees of Freedom
#'
#' This internal helper function computes the profile maximum likelihood estimates (MLEs)
#' of the regression coefficients (`beta`) and scale parameter (`sigma`) for a fixed
#' degrees of freedom (`nu`) under a Student-t error model.
#' It returns the estimated parameters along with the corresponding negative log-likelihood value
#' evaluated at these estimates.
#'
#' @param nu The fixed degrees of freedom for the t-distribution.
#' @param y The response variable.
#' @param x The design matrix including the intercept term.
#'
#' @return A list containing:
#' \describe{
#'   \item{beta_hat}{Profile MLE of `beta`.}
#'   \item{sigma_hat}{Profile MLE of `sigma`.}
#'   \item{loglik}{Negative log-likelihood value evaluated at the fitted parameters.}
#'   \item{convergence}{Convergence code (0 indicates success).}
#' }
#'
#' @keywords internal
#'
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

#' Profile MAP estimates of Sigma and Beta for Fixed Degrees of Freedom
#'
#' This function computes the MAP estimates
#' of regression coefficients (`beta`) and scale parameter (`sigma`) for a fixed
#' degrees of freedom value `nu` in a Student-t regression model.
#' The optimization is performed over the log-likelihood combined with a log-sigma
#' prior (i.e. log inverse prior) and log-beta prior (i.e. log(1) = 0) .
#'
#' This function serves as the first-stage optimization step in estimating the full
#' joint model parameters that include the prior on `nu` (i.e., log-likelihood +
#' log-sigma prior + log-beta prior + log-nu prior).
#'
#' @param nu The fixed degrees of freedom for the t-distribution.
#' @param y The response variable.
#' @param x The design matrix including the intercept term.
#'
#' @return A list containing:
#' \describe{
#'   \item{beta_hat}{Profile MLE of `beta`.}
#'   \item{sigma_hat}{Profile MLE of `sigma`.}
#'   \item{loglik}{Negative log-likelihood value evaluated at the fitted parameters.}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{beta_hat}{Estimated regression coefficients.}
#'   \item{sigma_hat}{Estimated scale parameter.}
#'   \item{loglik}{Maximized joint density (excluding `nu` prior).}
#'   \item{convergence}{Convergence code from \code{\link[stats]{optim}} (0 indicates success).}
#' }
#'
#' @keywords internal
#'
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
