# Check dimensions and inputs in MSAdata object

Ensures that data inputs are of proper dimension. Whenever possible,
default values are added to missing items.

## Usage

``` r
check_data(MSAdata, silent = FALSE)
```

## Arguments

- MSAdata:

  S4 object containing data inputs. See
  [MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)

- silent:

  Logical, whether or not to report default values to the console

## Value

An updated
[MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
object.

## See also

[MSAdata](https://blue-matter.github.io/multiSA/reference/MSAdata-class.md)
