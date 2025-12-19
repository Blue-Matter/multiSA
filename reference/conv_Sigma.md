# Calculate covariance matrix

Uses Cholesky factorization to generate a covariance matrix (or any
symmetric positive definite matrix).

## Usage

``` r
conv_Sigma(sigma, lower_diag)
```

## Arguments

- sigma:

  Numeric vector of marginal standard deviations (all greater than
  zeros)

- lower_diag:

  Numeric vector to populate the lower triangle of the correlation
  matrix. All real numbers. Length `sum(1:(length(sigma) - 1))`

## Value

Numeric

## Examples

``` r
set.seed(23)
n <- 5
sigma <- runif(n, 0, 2)
lower_diag <- runif(sum(1:(n-1)), -10, 10)
Sigma <- conv_Sigma(sigma, lower_diag)
Sigma/t(Sigma) # Is symmetric matrix? All ones
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    1    1    1    1
#> [2,]    1    1    1    1    1
#> [3,]    1    1    1    1    1
#> [4,]    1    1    1    1    1
#> [5,]    1    1    1    1    1
cov2cor(Sigma)
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  1.0000000 -0.8363414  0.6805096  0.7786780  0.6086421
#> [2,] -0.8363414  1.0000000 -0.1694648 -0.3245098 -0.3116680
#> [3,]  0.6805096 -0.1694648  1.0000000  0.9513182  0.6523565
#> [4,]  0.7786780 -0.3245098  0.9513182  1.0000000  0.7979121
#> [5,]  0.6086421 -0.3116680  0.6523565  0.7979121  1.0000000
```
