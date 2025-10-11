#' Simulate Linear Regression Data with t-Distributed Errors
#'
#' This function generates synthetic data from a linear regression model
#' where the error terms follow a Student-t distribution with specified
#' degrees of freedom. It is useful for studying the performance of robust
#' regression methods under heavy-tailed or outlier-prone error distributions.
#
#' The model is defined as:
#' \deqn{y_i = x_i^\top \beta + \varepsilon_i, \quad \varepsilon_i \sim t_\nu(0, \sigma^2)}
#'
#' @param n Number of observations to simulate. Default is 300.
#' @param p Number of predictors (including intercept). Default is 1.
#' @param beta Numeric vector of regression coefficients of length \code{p}.
#'   The first element corresponds to the intercept.
#' @param sigma Scale parameter for the t-distributed errors. Default is 1.
#' @param nu Degrees of freedom for the t errors. Default is 2.
#' @param seed Optional seed
#'
#' @return A list with the following components:
#' \describe{
#'   \item{y}{A numeric vector of simulated response values.}
#'   \item{x}{A numeric matrix of predictors (including intercept) from standard normal.}
#' }
#'
#' @examples
#' # Simulate data with heavy-tailed errors
#' sim <- simulate_t_error_data(
#'   n = 300, p = 2, beta = c(1, 2), sigma = 1, nu = 3, seed = 123
#' )
#' head(sim$x)
#' head(sim$y)
#'
#' @export
simulate_t_error_data <- function(n = 300, p = 1, beta = rep(0,1), sigma = 1, nu = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(beta) != p) stop("length(beta) must equal p")
  ncov <- p - 1
  x <- matrix(0, nrow = n, ncol = ncov)
  if (ncov > 0) for (i in 1:ncov) x[,i] <- rnorm(n)
  x <- cbind(1, x)
  mu <- x %*% beta
  errors <- sigma * rt(n = n, df = nu)
  y <- as.numeric(mu + errors)
  list(y = y, x = x)
}
