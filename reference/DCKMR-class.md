# DCKMR S4 object

Store lists of data frames for parent offspring pairs and half-sibling
pairs

## Slots inherited from DCKMR

- `POP_s`:

  A list by stock of data frames for parent-offspring pairs. Each row in
  the data frame corresponds to a "sampling unit" defined by the
  columns:

  |     |                                |
  |-----|--------------------------------|
  | `a` | Capture year of parent         |
  | `t` | Age at capture of parent       |
  | `y` | Birth year of offspring        |
  | `n` | Number of pairwise comparisons |
  | `m` | Number of POPs                 |

- `HSP_s`:

  A list by stock of data frames for half-sibling pairs. Each row in the
  data frame corresponds to a "sampling unit" defined by the columns:

  |      |                                |
  |------|--------------------------------|
  | `yi` | Birth year of older sibling    |
  | `yj` | Birth year of younger sibling  |
  | `n`  | Number of pairwise comparisons |
  | `m`  | Number of HSPs                 |

- `CKMR_like`:

  Character, likelihood for the POP and HSP sampling units. See `type`
  argument in
  [`like_CKMR()`](https://blue-matter.github.io/multiSA/reference/like_CKMR.md)
  for options.
