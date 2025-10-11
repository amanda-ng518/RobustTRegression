#' Estimate Regression Coefficients \eqn{\beta} Using Various Methods
#'
#' This function estimates regression coefficients (\eqn{\beta}) and the scale
#' parameter (\eqn{\sigma}) using one of several estimation approaches:
#'
#' \itemize{
#'   \item Ordinary Least Squares: \code{"OLS"}
#'   \item Huber Regression (robust): \code{"Huber"}
#'   \item Fixed degrees of freedom Student-\emph{t} regression: a non-negative integer \code{nu}
#'   \item Estimated degrees of freedom using one of the following:
#'     \code{"profile"}, \code{"adj_profile"}, \code{"IJ"}, \code{"Marginalized IJ"}, \code{"Marginalized Fisher"}
#' }
#'
#' Internally, this function dispatches to the appropriate estimation routine based on the specified method.
#'
#' @param y Numeric vector. The response variable.
#' @param x Numeric matrix. The design matrix (must include an intercept column).
#' @param method Character or numeric. Specifies the estimation method:
#'   one of \code{"OLS"}, \code{"Huber"}, \code{"profile"}, \code{"adj_profile"},
#'   \code{"IJ"}, \code{"Marginalized IJ"}, \code{"Marginalized Fisher"},
#'   or a positive integer (interpreted as a fixed \code{nu}).
#' @param omega_init Numeric. Initial value for the \code{omega = 1/nu} parameter (default is 0.5).
#' @param beta_init Optional numeric vector. Initial values for \code{beta}.
#'   If not specified, ordinary least squares (OLS) estimates are used as defaults.
#' @param sigma_init Optional numeric value. Initial value for \code{sigma}.
#'   If not specified, the residual standard deviation from OLS is used as the default.
#' @param control Optional list of control parameters passed to \code{\link[stats]{optim}}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{method}{Character. The estimation method used.}
#'   \item{nu_hat}{Estimated or fixed degrees of freedom \eqn{\nu}. \code{NA} for OLS and Huber.}
#'   \item{omega_hat}{Estimated or fixed value of \eqn{\omega = 1/\nu}. \code{NA} for OLS and Huber.}
#'   \item{beta_hat}{Estimated regression coefficients.}
#'   \item{sigma_hat}{Estimated scale parameter (residual standard deviation).}
#'   \item{success_nu}{Convergence code for \eqn{\nu} estimation. \code{0} indicates successful convergence; \code{NA} for OLS, Huber, and fixed \code{nu}.}
#'   \item{success_beta}{Convergence code for \code{beta}/\code{sigma} optimization. \code{0} indicates success; \code{NA} for OLS and Huber.}
#' }
#'
#' @importFrom MASS rlm psi.huber
#' @export
estimate_beta <- function(y, x,
                          method = c("OLS", "Huber", "Profile", "Adj profile", "IJ",
                                     "Marginalized IJ", "Marginalized Fisher"),
                          omega_init = 0.5,
                          beta_init = NULL,
                          sigma_init = NULL,
                          control = list(maxit = 10000)) {
  method <- as.character(method)
  n <- length(y)
  p <- ncol(x)

  # ---- CASE 1: OLS ----
  if (method == "OLS") {
    lmfit <- lm(y ~ x - 1)
    beta_hat <- coef(lmfit)
    sigma_hat <- sd(residuals(lmfit))
    return(list(
      method = method,
      nu_hat = NA,
      omega_hat = NA,
      beta_hat = beta_hat,
      sigma_hat = sigma_hat,
      success_nu = NA,
      success_beta = NA
    ))
  }

  # ---- CASE 2: Huber regression ----
  if (method == "Huber") {
    fit_lm <- lm(y ~ x - 1)
    RB <- rlm(y ~ x[, -1], psi = psi.huber, k = 1.345, maxit = 100)
    beta_huber <- rlmDD(y, x[, -1], fit_lm$coef[-2], RB$coef,
                        method = "Huber", plot = "N")$esti[["coefficients"]]
    sigma_hat <- sd(y - as.vector(x %*% beta_huber))
    return(list(
      method = method,
      nu_hat = NA,
      omega_hat = NA,
      beta_hat = beta_huber,
      sigma_hat = sigma_hat,
      success_nu = NA,
      success_beta = NA
    ))
  }

  # ---- CASE 3: Fixed nu or nu-estimation method ----

  if (is.numeric(method)) method <- as.character(method)
  valid_methods <- c("Profile", "Adj profile", "IJ", "Marginalized IJ", "Marginalized Fisher")

  # 3a. Estimate nu if method is one of the 5 estimation methods
  if (method %in% valid_methods) {
    est <- switch(method,
                  "Profile" = estimate_nu_profile(y, x, omega_init),
                  "Adj profile" = estimate_nu_adj_profile(y, x, omega_init),
                  "IJ" = estimate_nu_IJ(y, x, omega_init),
                  "Marginalized IJ" = estimate_nu_mar_IJ(y, x, omega_init),
                  "Marginalized Fisher" = estimate_nu_nu_block(y, x, omega_init)
    )
    omega_hat <- est$omega
    nu_hat <- ifelse(omega_hat > 0, 1 / omega_hat, Inf)
    success_nu <- est$convergence
  } else {
    # 3b. Fixed nu
    if (suppressWarnings(!is.na(as.numeric(method))) &&
               as.numeric(method) > 0) {
      nu_hat <- as.numeric(method)
      omega_hat <- ifelse(nu_hat == 0, Inf, 1 / nu_hat)
    } else {
      stop("Invalid method: must be 'OLS', 'Huber', one of the 5 nu estimation methods, or a non-negative integer.")
    }
    success_nu <- NA
  }

  # 3c. Optimize beta and sigma given nu
  fit_res <- fit_profile_mle_fixed_nu(y, x, nu_hat, par_init = NULL, control = control)

  list(
    method = method,
    nu_hat = nu_hat,
    omega_hat = omega_hat,
    beta_hat = fit_res$beta,
    sigma_hat = fit_res$sigma,
    success_nu = success_nu,
    success_beta = fit_res$convergence
  )
}
