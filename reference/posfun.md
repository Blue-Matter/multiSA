# Quadratic penalty function

Taped penalty function if `x < eps`

## Usage

``` r
posfun(x, eps)
```

## Arguments

- x:

  Numeric, the parameter

- eps:

  Numeric, the threshold below which a penalty will be applied

## Value

The penalty value is

\$\$ \textrm{penalty} = \begin{cases} 0.1 (x - \varepsilon)^2 & x \le
\varepsilon\\ 0 & x \> \varepsilon \end{cases} \$\$

Numeric
