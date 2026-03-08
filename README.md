# RobustTRegression

**RobustTRegression** is an R package for performing robust regression using the t-distribution, as discussed in [Ng et al. (2026)](https://arxiv.org/abs/2603.00269). It provides tools to fit regression models that are resistant to outliers and heavy-tailed errors, improving the reliability of parameter estimates in such scenarios. For a more complete demonstration to reproduce the numerical results in our paper, please refer to our full GitHub repo: [Aspect of Robust Regression Analysis](https://github.com/amanda-ng518/Aspects-of-Robust-Regression-Analysis).

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

### Data simulations

To generate synthetic data from a linear regression model where the error terms follow a Student's t-distribution with a specified degrees of freedom (`nu`), use `simulate_t_error_data()`. Covariates x are generated from standard normal distribution. Optionally, you may set a random `seed` for reproducibility. By default, `sigma = 1, nu = 2`. 

Example:
```r
# Load library
> library(RobustTRegression)

# Simulate a set of data with 4 rows, 5 columns with true beta all set as 0 and t(2) error
> sim_t_data <- simulate_t_error_data(n = 4, p = 5, beta = rep(0,5), sigma = 1, nu = 2, seed = 123)

> sim_t_data$x
     [,1]        [,2]       [,3]       [,4]       [,5]
[1,]    1 -0.56047565  0.1292877 -0.6868529  0.4007715
[2,]    1 -0.23017749  1.7150650 -0.4456620  0.1106827
[3,]    1  1.55870831  0.4609162  1.2240818 -0.5558411
[4,]    1  0.07050839 -1.2650612  0.3598138  1.7869131

> sim_t_data$y
[1]  0.9913147 -2.0389531 -1.0818755  0.7809687
```

To generate synthetic data from a linear regression model where the error terms follow a normal distribution with a specified mean (`mean`) and standard error (`sigma`), use `simulate_n_error_data()`. Covariates x are generated from standard normal distribution. Optionally, you may set a random `seed` for reproducibility. By default, `mean = 0, sigma = 1`.

Example:
```r
# Simulate a set of data with 4 rows, 5 columns with true beta all set as 0 and standard normal error
sim_n_data <- simulate_n_error_data(n = 4, p = 5, beta = rep(0,5), mean = 0, sigma = 1, seed = 123)
x <- sim_n_data$x
y <- sim_n_data$y
```

To generate synthetic data from a linear regression model with normal errors (specify `mean` and `sigma`) and contamination, use `simulate_contaminated_data()`. Covariates x are generated from standard normal distribution. The contamination is controlled by `contam_type` and `contam_prob`. 

The available contamination types in `contam_type` are:

- "N_0_9" : $N(0,9)$ error
- "t_2" : $t(2)$ error
- "chisq" : $\chi^2(4)$ - 4 error
- "twopt" : two-point contamination in which any response value may become $-5$ or 5 with a probability $\lambda /2$ each, where $\lambda$ is the contamination probability.
- `NULL` (also set `contam_prob` = 0): ignore contamination specification

`contam_prob` should be a number between 0 to 1, controlling for the probability of contamination in the data. Remember to choose a non-null option in `contam_type` if `contam_prob` is set as non-zero. Optionally, you may a set random `seed` for reproducibility. By default, `mean = 0, sigma = 1`.

Example:
```r
# Simulate a set of data with 20 rows, 5 columns with true beta all set as 0 and t(2) error
contam_sim_data <- simulate_contaminated_data(n = 20, p = 5, beta = rep(0, 5), mean = 0,
                                       sigma = 1, contam_type = "twopt", contam_prob = 0.2,
                                       seed = 123)
x <- contam_sim_data$x
y <- contam_sim_data$y
```

### $\nu$-estimation

$\nu$ estimation objective functions:

- `estimate_nu_profile()`: Profile likelihood  
- `estimate_nu_adj_profile()`: Adjusted profile likelihood
- `estimate_nu_IJ()`: Full Bayes
- `estimate_nu_nu_block()`: Pseudo Bayes

For each of the 4 $\nu$ estimation approaches, we provide functions to 
compute the 4 objective functions estimator for the degrees of freedom
parameter $\nu$ via optimization over $\omega = 1 / \nu$. Users need to
specify `x`, `y`, and the initial value in $\omega$ scale for BFGS optimization. By default, 
it is set as 0.5 (i.e., $\nu = 2$).

They output a list of: 

- `omega` : optimal $\omega$
- `nu` : optimal $\nu$ (essentially 1/ optimal $\omega$)
- `convergence` : code to check BFGS convergence (0 if success)
  
Example:
```r
# Estimate nu with profile likelihood approach
# Using data in `contam_sim_data` generated above
> nu_estimation <- estimate_nu_profile(y, x, omega_init = 1/2)
> nu_estimation$nu # Optimal nu
[1] 0.1102719
> nu_estimation$convergence # Check convergence = 0 (if success)
[1] 0
```

We also provide a function to conduct estimations on all 4 approaches at once. Similarly, users need to specify `x`, `y`, and the initial value in $\omega$ scale for BFGS optimization. It outputs a list of 4 sub-lists (`profile`, `adj_profile`, `IJ`, `nu_block`), each sub-list containing the results from a single estimation method (corresponding to profile, adjusted profile, full bayes, and pseudo bayes).

Example:
```r
# Using data in `contam_sim_data` generated above
> run_all_estimators(y, x, omega_init = 1/2)
$profile
$profile$omega
[1] 9.068491

$profile$nu
[1] 0.1102719

$profile$convergence
[1] 0

$profile$value
[1] 26.64138


$adj_profile
$adj_profile$omega
[1] 0.7149223

$adj_profile$nu
[1] 1.398753

$adj_profile$convergence
[1] 0

$adj_profile$value
[1] 48.60885


$IJ
$IJ$omega
[1] 7.822529

$IJ$nu
[1] 0.1278359

$IJ$convergence
[1] 0

$IJ$value
[1] 16.26669


$nu_block
$nu_block$omega
[1] 9.292356

$nu_block$nu
[1] 0.1076153

$nu_block$convergence
[1] 0

$nu_block$value
[1] 28.35427
```

### $\beta$-estimation

We provide function `estimate_beta()` to first compute a $\nu$ estimator and then proceed to 
optimize the t-likelihood to obtain $\beta$ estimators. Aside from `x`, `y` and starting point `omega_init`, you will need to specify the `method` to conduct $\nu$ estimation. The available methods include:

- "OLS" : $\nu$ will not be estimated
- "Huber" : $\nu$ will not be estimated
- "Profile"
- "Adj profile"
- "Full Bayes"
- "Pseudo Bayes"
- any positive integer (interpreted as a fixed nu) : $\hat \nu$ will be set as this fixed integer

It outputs a list of : 

- `method` : method used to estimate `nu`, same as your input `method`
- `nu_hat` : optimal $\nu$ (essentially 1/ optimal $\omega$)
- `beta_hat` : optimal $\beta$
- `sigma_hat` : optimal $\sigma$
- `success_nu` : Convergence code for $\nu$-estimation (0 if success)
- `success_beta` : Convergence code for $\beta$-estimation (0 if success)

Example:
```r
# Estimate beta with adjusted profile likelihood approach
# Using data in `contam_sim_data` generated above
> beta_estimation <- estimate_beta(y, x, method = "Adj profile", omega_init = 0.5) # By default we use OLS estimates as initial guess for beta and sigma
> beta_estimation$beta_hat # Optimal beta
         x1          x2          x3          x4          x5 
 0.24335350  0.23649201  0.23051357 -0.03220889 -0.09397384 
> beta_estimation$nu_hat # Optimal nu
[1] 1.398753
> beta_estimation$success_beta # Check convergence = 0 (if success)
[1] 0
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

# Load library
library(RobustTRegression)
```
