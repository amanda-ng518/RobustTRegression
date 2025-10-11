# RobustTRegression

**RobustTRegression** is an R package for performing robust regression using the t-distribution. It provides tools to fit regression models that are resistant to outliers and heavy-tailed errors, improving the reliability of parameter estimates in such scenarios.

## Features

- Robust regression modeling using t-distributed errors
- Parameter estimation resistant to outliers
- Easy-to-use interface similar to standard regression functions
- Diagnostic tools to assess model fit

## Installation

You can install the package directly from GitHub using `remotes`:

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install RobustTRegression from GitHub
remotes::install_github("amanda-ng518/RobustTRegression")
