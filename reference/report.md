# Generate markdown reports

Generate a markdown report of model fits and estimates.

## Usage

``` r
report(object, ...)

# S4 method for class 'MSAassess'
report(
  object,
  name,
  filename = "MSA",
  dir = tempdir(),
  open_file = TRUE,
  render_args = list(),
  ...
)
```

## Arguments

- object:

  An object from MSA.

- ...:

  Additional arguments to render reports.

- name:

  Optional character string for the model name to include in the report,
  e.g., model run number. Default uses `substitute(object)`

- filename:

  Character string for the name of the markdown and HTML files.

- dir:

  The directory in which the markdown and HTML files will be saved.

- open_file:

  Logical, whether the HTML document is opened after it is rendered.

- render_args:

  List of arguments to pass to
  [`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html).

## Value

`report` generic for MSAassess invisibly returns the output of
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html):
character of the path of the rendered HTML markdown report.
