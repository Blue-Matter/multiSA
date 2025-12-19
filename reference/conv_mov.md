# Calculate movement matrix for all age classes

Movement matrices are calculated for all age classes from a base matrix
and a gravity model formulation (Carruthers et al. 2016).

## Usage

``` r
conv_mov(x, g, v, na = dim(x)[1], nr = dim(x)[2], aref = ceiling(0.5 * na))
```

## Arguments

- x:

  Base log-movement parameters. See details. Array `[a, r, r]`

- g:

  Gravity model attractivity term. Tendency to move to region `r`.
  Matrix `[a, r]`

- v:

  Gravity model viscosity term. Tendency to stay in same region. Vector
  by `a`

- na:

  Integer, number of ages

- nr:

  Integer, number of regions

- aref:

  Integer, reference age class

## Value

Movement array `[a, r, r]`

## Details

Rows index region of origin and columns denote region of destination.

In log space, the movement matrix \\m_a\\ for age class \\a\\ from
region \\r\\ to \\r'\\ is the sum of base matrix \\x\\ and gravity
matrix \\G\\: \$\$m\_{a,r,r'} = x\_{a,r,r'} + G\_{a,r,r'}\$\$

To essentially exclude movement from \\r\\ to \\r'\\, set \\x\_{a,r,r'}
= -1000\\.

Gravity matrix \\G\\ includes an attractivity term \\g\\ and viscosity
term \\v\\:

\$\$G\_{a,r,r'} = \begin{cases} g'\_{a,r'} + v_a \quad & r = r'\\
g'\_{a,r'} \quad & \textrm{otherwise} \end{cases} \$\$

Vector \\g'\\ are offset terms relative to the value for the reference
age class: \$\$g'\_{a,r'} = \begin{cases} g\_{a,r} \quad & a =
a\_{ref}\\ g\_{a,r} + g\_{a=aref,r} \quad & \textrm{otherwise}
\end{cases} \$\$

The movement matrix in normal space is obtained by the softmax
transformation: \$\$M\_{a,r,r'} =
\dfrac{\exp(m\_{a,r,r'})}{\sum\_{r'}\exp(m\_{a,r,r'})}\$\$

If \\x\\ and \\v\\ are zero, then the movement matrix simply distributes
the total stock abundance into the various regions as specified in
\\g'\\.

## References

Carruthers, T.R., et al. 2015. Modelling age-dependent movement: an
application to red and gag groupers in the Gulf of Mexico. CJFAS 72:
1159-1176.
[doi:10.1139/cjfas-2014-0471](https://doi.org/10.1139/cjfas-2014-0471)
