# RobustTRegression

**RobustTRegression** is an R package for performing robust regression using the t-distribution. It provides tools to fit regression models that are resistant to outliers and heavy-tailed errors, improving the reliability of parameter estimates in such scenarios.

## Overview

We utilize two-stage procedure under the student’s-t regression framework
as follows: we first estimate the degrees of freedom parameter ($\nu$), and then
use this estimate to infer the regression coefficients ($\beta$). We provide several 
frequentist and bayesian methods, including:

- Profile likelihood
- Adjusted profile likelihood
- Full Bayes
- Pseudo Bayes

## Functions

### Simulate data
We provide 
Generate synthetic data from a linear regression model where the error terms
follow a Student's t-distribution with a specified degrees of freedom.

Example:
```r
simulate_t_error_data(n = 300, p = 1, beta = rep(0,1), sigma = 1, nu = 2, seed = NULL)
```

### Estimate $\nu$

$\nu$ estimation objective functions:

- `estimate_nu_profile()`: Profile likelihood  
- `estimate_nu_adj_profile()`: Adjusted profile likelihood
- `estimate_nu_IJ()`: Full Bayes
- `estimate_nu_nu_block()`: Pseudo Bayes

For each of the 4 $\nu$ estimation approaches, we provide functions to 
compute the 4 objective functions estimator for the degrees of freedom
parameter \eqn{\nu} via optimization over \eqn{\omega = 1 / \nu}. Users need to
specify `x`, `y`, and the initial value in $\omega$ scale for BFGS optimization. By default, 
it is set as 0.5 (i.e., $\nu = 2$).

They output a list of: 

- `omega`: optimal $\omega$
- `nu`: optimal $\nu$ (essentially 1/ optimal $\omega$)
- `convergence`: code to check BFGS convergence (0 if success)
  
Example:
```r
# Estimate nu with profile likelihood approach
nu_estimation <- estimate_nu_profile(y, x, omega_init = 1/2)
est_nu <- nu_estimation$nu # Optimal nu
nu_estimation$convergence # Check convergence = 0 (if success)
```

We also provide a function to conduct estimations on all 4 approaches at once. Similarly, users need to specify `x`, `y`, and the initial value in $\omega$ scale for BFGS optimization. It outputs a list of 4 sub-lists, each sub-list containing the results from a single estimation method.

Example:
```r
run_all_estimators(y, x, omega_init = 1/2)
```

### Estimate $\beta$

We provide function `estimate_beta()` to first compute a $\nu$ estimator and then proceed to 
optimize the t-likelihood to obtain $\beta$ estimators. Aside from `x`, `y` and starting point `omega_init`, you will need to specify the `method` to conduct $\nu$ estimation. The available methods include:

- "OLS": $\nu$ will not be estimated
- "Huber": $\nu$ will not be estimated
- "Profile"
- "Adj profile"
- "Full Bayes"
- "Pseudo Bayes"
- any positive integer (interpreted as a fixed nu): $\hat \nu$ will be set as this fixed integer

It outputs a list of : 

- `method`: method used to estimate `nu`, same as your input `method`
- `nu_hat`: optimal $\nu$ (essentially 1/ optimal $\omega$)
- `beta_hat`: optimal $\beta$
- `sigma_hat`: optimal $\sigma$
- `success_nu`: Convergence code for $\nu$-estimation (0 if success)
- `success_beta`: Convergence code for $\beta$-estimation (0 if success)

Example:
```r
# Estimate beta with adjusted profile likelihood approach
beta_estimation <- estimate_beta(y, x, method = "Adj profile", omega_init = 0.5) # By default we use OLS estimates as initial guess for beta and sigma
est_beta <- beta_estimation$beta_hat # Optimal beta
nu_hat <- beta_estimation$nu_hat # Optimal nu
beta_estimation$success_beta # Check convergence = 0 (if success)
```


## Installation

You can install the R package directly from GitHub using `remotes`:

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install RobustTRegression from GitHub
remotes::install_github("amanda-ng518/RobustTRegression")
```
