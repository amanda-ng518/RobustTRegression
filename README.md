# RobustTRegression

**RobustTRegression** is an R package for performing robust regression using the t-distribution. It provides tools to fit regression models that are resistant to outliers and heavy-tailed errors, improving the reliability of parameter estimates in such scenarios.

## Overview

We utilize two-stage procedure under the student’s-t regression framework
as follows: we first estimate the degrees of freedom parameter ($\nu$), and then
use this estimate to infer the regression coefficients ($\beta$). We provide several 
frequentist and bayesian methods, including:

- Profile likelihood
- Adjusted profile likelihood
- Independent Jeffrey's 
- Marginalized Jeffrey's
- Marginalized Fisher

## Functions

### Simulate data
Generate synthetic data from a linear regression model where the error terms
follow a Student's t-distribution with a specified degrees of freedom.

Example:
```r
simulate_t_error_data(n = 300, p = 1, beta = rep(0,1), sigma = 1, nu = 2, seed = NULL)
```

### Estimate $\nu$

$\nu$ estimation objective functions:

- Profile likelihood
- Adjusted profile likelihood
- Likelihood with sigma inverse prior, constant $\beta$ prior and $\nu$ Independent Jeffrey's prior
- Profile likelihood with $\nu$ Jeffrey's prior
- Profile likelihood with $\nu$ Fisher prior ($\nu$ block in Fisher information)

For each of the 5 $\nu$ estimation approaches, we provide functions to 
compute the 5 objective functions estimator for the degrees of freedom
parameter \eqn{\nu} via optimization over \eqn{\omega = 1 / \nu}. Users need to
specify the initial value in $\omega$ scale for BFGS optimization. By default, 
it is set as 0.5 (i.e., $\nu = 2$).

Example:
````r
estimate_nu_profilefunction(y, x, omega_init = 1/2)
```

We also provide a function to conduct estimations on all 5 approaches all at once.
Users need to specify the initial value in $\omega$ scale for BFGS optimization.

Example:
````r
run_all_estimators(y, x, omega_init = 1/2)
```

### Estimate $\beta$

We provide functions to first compute a $\nu$ estimator and then proceed to 
optimize the t-likelihood to obtain $\beta$ estimators. Users can choose one 
of the following $\nu$-estimation approaches:

- Profile likelihood
- Adjusted profile likelihood
- Independent Jeffrey's 
- Marginalized Jeffrey's 
- Marginalized Fisher 

Users need to specify the initial value in $\omega$ scale for BFGS optimization.
By default, it is set as 0.5 (i.e., $\nu = 2$).

Additionally, users can also opt for non $\nu$-estimation approaches, which
prodcue $\beta$ estimators directly:

- Fixed $\nu$ (specify `method` as a positive integer)
- OLS
- Huber

Example:
````r
estimate_beta(y, x, method = "Profile", omega_init = 0.5)
```


## Installation

You can install the package directly from GitHub using `remotes`:

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
```

```r
# Install RobustTRegression from GitHub
remotes::install_github("amanda-ng518/RobustTRegression")
```
