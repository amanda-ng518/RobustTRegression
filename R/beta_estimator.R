# ============================================================
# Beta estimation function
# ============================================================
#' Estimate regression coefficients beta under various methods
#'
#' This function estimates regression coefficients (\eqn{\beta}) and scale (\eqn{\sigma})
#' using one of the following approaches:
#' - Ordinary least squares ("OLS")
#' - Huber regression ("Huber")
#' - Fixed degrees of freedom t-regression (method = non-negative integer)
#' - Estimated degrees of freedom using one of five methods:
#'   "Profile", "Adj profile", "IJ", "Marginalized IJ", or "Marginalized Fisher"
#'
#' @importFrom MASS rlm psi.huber
#'
#' @param y Response vector
#' @param x Design matrix (including intercept)
#' @param method Character or numeric; one of:
#'   "OLS", "Huber", "profile", "adj_profile", "IJ", "Marginalized IJ", "Marginalized Fisher",
#'   or a non-negative integer for fixed nu approach.
#' @param omega_init Initial value for omega (default = 0.5)
#' @param beta_init Optional initial beta values
#' @param sigma_init Optional initial sigma
#' @param control List of control parameters for optim
#'
#' @return A list with components:
#' \describe{
#'   \item{method}{The method used}
#'   \item{nu_hat}{Estimated or fixed degrees of freedom (NA for OLS/Huber)}
#'   \item{omega_hat}{Transformed 1/nu parameter (NA for OLS/Huber)}
#'   \item{beta_hat}{Estimated coefficients}
#'   \item{sigma_hat}{Estimated sigma (residual SD)}
#'   \item{success_nu}{Convergence code for nu estimation. A value of
#'     0 indicates successful convergence. (NA for fixed nu/OLS/Huber)}
#'   \item{success_beta}{Convergence code for beta optimization. A value of
#'     0 indicates successful convergence. (NA for OLS/Huber)}
#' }
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
               as.numeric(method) >= 0) {
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
