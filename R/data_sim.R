# ============================================================
# Data Simulation Functions
# ============================================================

#' Simulate Linear Regression Data with t-Distributed Errors
#'
#' Generate synthetic data from a linear regression model where the error terms
#' follow a Student's t-distribution with a specified degrees of freedom.
#' This is useful for studying the performance of robust regression methods
#' under heavy-tailed or outlier-prone error distributions.
#'
#' The model is defined as:
#' \deqn{y_i = x_i^\top \beta + \varepsilon_i, \quad \varepsilon_i \sim t_\nu(0, \sigma^2)}
#'
#' @param n Integer. Number of observations to simulate. Default is 300.
#' @param p Integer. Number of predictors (including the intercept). Default is 1.
#' @param beta Numeric vector of regression coefficients of length \code{p}.
#' The first element corresponds to the intercept.
#' @param sigma Numeric. Scale parameter for the t-distributed errors. Default is 1.
#' @param nu Numeric. Degrees of freedom for the t-distribution. Default is 2.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{y}{A numeric vector of simulated response values.}
#'   \item{x}{A numeric matrix of predictor values (including intercept),
#'   drawn from a standard normal distribution.}
#' }
#'
#' @examples
#' set.seed(123)
#' data <- simulate_t_error_data(n = 100, p = 3, beta = c(1, 0.5, -0.3), sigma = 1, nu = 3)
#' str(data)
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


#' Simulate Linear Regression Data with Normal Errors
#'
#' Generate synthetic data from a linear regression model where the error terms
#' follow a normal distribution.
#' This is useful for studying the performance of robust regression methods
#' under heavy-tailed or outlier-prone error distributions.
#'
#' The model is defined as:
#' \deqn{y_i = x_i^\top \beta + \varepsilon_i, \quad \varepsilon_i \sim N(\mu, \sigma^2)}
#'
#' @param n Integer. Number of observations to simulate. Default is 300.
#' @param p Integer. Number of predictors (including the intercept). Default is 1.
#' @param beta Numeric vector of regression coefficients of length \code{p}.
#' The first element corresponds to the intercept.
#' @param mean Numeric. Mean parameter for the normal errors. Default is 0.
#' @param sigma Numeric. Scale parameter for the normal errors. Default is 1.
#' @param seed Optional integer. Set random seed for reproducibility.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{y}{A numeric vector of simulated response values.}
#'   \item{x}{A numeric matrix of predictor values (including intercept),
#'   drawn from a standard normal distribution.}
#' }
#'
#' @examples
#' set.seed(123)
#' data <- simulate_n_error_data(n = 100, p = 3, beta = c(1, 0.5, -0.3), mean = 0, sigma = 1)
#' str(data)
#'
#' @export
simulate_n_error_data <- function(n = 300, p = 1, beta = rep(0,1), mean = 0, sigma = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(beta) != p) stop("length(beta) must equal p")
  ncov <- p - 1
  x <- matrix(0, nrow = n, ncol = ncov)
  if (ncov > 0) for (i in 1:ncov) x[,i] <- rnorm(n)
  x <- cbind(1, x)
  mu <- x %*% beta
  errors <- rnorm(n = n, mean = mean, sd = sigma)
  y <- as.numeric(mu + errors)
  list(y = y, x = x)
}

#' Simulate Linear Regression Data with Contaminated Errors
#'
#' Generate synthetic data from a linear regression model where the error terms
#' follow a normal distribution, with prespecified contamination.
#' This is useful for studying the performance of robust regression methods
#' under heavy-tailed or outlier-prone error distributions.
#'
#' The model is defined as:
#' \deqn{y_i = x_i^\top \beta + \varepsilon_i, \quad \varepsilon_i \sim N(\mu, \sigma^2)}
#'
#' @param n Integer. Number of observations to simulate. Default is 300.
#' @param p Integer. Number of predictors (including the intercept). Default is 1.
#' @param beta Numeric vector of regression coefficients of length \code{p}.
#' The first element corresponds to the intercept.
#' @param mean Numeric. Mean parameter for the normal errors. Default is 0.
#' @param sigma Numeric. Scale parameter for the normal errors. Default is 1.
#' @param contam_type Character. Type of contaminated error. Default is NULL.
#' Available options include: "N_0_9", "t_2", "chisq", "twopt".
#' @param contam_prob Numeric. Proportion of data with contaminated error. Default is 0.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{y}{A numeric vector of simulated response values.}
#'   \item{x}{A numeric matrix of predictor values (including intercept),
#'   drawn from a standard normal distribution.}
#' }
#'
#' @examples
#' set.seed(123)
#' data <- simulate_contaminated_data(n = 100, p = 3, beta = c(1, 0.5, -0.3), mean = 0, sigma = 1, contam_type = "twopt", contam_prob = 0.1)
#' str(data)
#'
#' @export
simulate_contaminated_data <- function(n = 300, p = 1, beta = rep(0, 1), mean = 0,
                                       sigma = 1, contam_type = NULL, contam_prob = 0,
                                       seed = NULL){
  if (!is.null(seed)) set.seed(seed)
  if (length(beta) != p) stop("length(beta) must equal p")
  if (is.null(contam_type) & contam_prob != 0) stop("Please specify contam_type")
  if (contam_prob < 0 | contam_prob > 1) stop("contam_prob must be between 0 and 1")

  # Design matrix
  ncov <- p - 1
  x <- matrix(0, nrow = n, ncol = ncov)
  if (ncov > 0) for (i in 1:ncov) x[, i] <- rnorm(n)
  x <- cbind(1, x)

  mu <- as.numeric(x %*% beta)

  # Initialize errors
  errors <- numeric(n)

  # Number of contaminated points
  m <- floor(n * contam_prob)

  # Contamination types
  # No contamination
  if (is.null(contam_type)){
    errors <- rnorm(n, mean = mean, sd = sigma)
  }
  # N(0,9)
  else if (contam_type == "N_0_9") {
    errors[1:m] <- rnorm(m, mean = 0, sd = 3)
    errors[(m+1):n] <- rnorm(n - m, mean = mean, sd = sigma)

  # t(2)
  } else if (contam_type == "t_2") {
    errors[1:m] <- sigma * rt(m, df = 2)
    errors[(m+1):n] <- rnorm(n - m, mean = mean, sd = sigma)

  # chisq(4) - 4
  } else if (contam_type == "chisq") {
    errors[1:m] <- rchisq(m, df = 4) - 4
    errors[(m+1):n] <- rnorm(n - m, mean = mean, sd = sigma)

  # two-point contamination
  } else if (contam_type == "twopt") {
    half <- floor(m / 2)
    y <- mu
    y[1:half] <- 5
    y[(half+1):m] <- -5
    y[(m+1):n] <- mu[(m+1):n] + rnorm(n - m, mean = mean, sd = sigma)
    return(list(y = y, x = x))

  } else {
    stop("Please specify correct contiminated error, or leave it as NULL")
  }

  y <- mu + errors
  list(y = y, x = x)
}
