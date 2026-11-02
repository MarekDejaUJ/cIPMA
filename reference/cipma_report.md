# Write a Full cIPMA Report

Writes the full summary of a cIPMA, and of its diagnostics when given,
as a Markdown report.

## Usage

``` r
cipma_report(x, file = NULL, diagnostics = NULL)
```

## Arguments

- x:

  A cipma object.

- file:

  Optional path of the Markdown file; the report is returned only when
  `file` is `NULL`.

- diagnostics:

  Optional result of
  [`cipma_diagnostics`](https://marekdejauj.github.io/cIPMA/reference/cipma_diagnostics.md).

## Value

The lines of the report, invisibly.

## Examples

``` r
# \donttest{
library(seminr)
mm <- constructs(
  composite("Image", multi_items("IMAG", 1:5)),
  composite("Value", multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3))
)
sm <- relationships(paths(from = c("Image", "Value"), to = "Satisfaction"))
pls <- estimate_pls(mobi, mm, sm)
#> Generating the seminr model
#> All 250 observations are valid.
boot <- bootstrap_model(pls, nboot = 50, cores = 1, seed = 1)
#> Bootstrapping model using seminr...
#> SEMinR Model successfully bootstrapped
scales <- setNames(rep(list(c(1, 10)), ncol(mobi)), names(mobi))
res <- cipma(boot, "Satisfaction", mobi, scales, target_level = 80, B = 99, seed = 1)
report <- cipma_report(res, file = file.path(tempdir(), "cipma-report.md"))
# }
```
