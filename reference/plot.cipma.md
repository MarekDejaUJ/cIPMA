# Plot cIPMA Results

Draws the combined map, the tolerance path of the effect sizes or the
tolerance path of the required levels, in grey scale. On the map, white
bubbles mark the predictors necessary at `pi = 1`, with a size given by
the share of cases below the required level at the chosen tolerance, and
black dots mark the other predictors. `scale = "zoom"` fits the axes to
the results; `scale = "full"` shows performance on the whole 0-100 scale
and importance from 0 to the next tenth above the largest total effect.

## Usage

``` r
# S3 method for class 'cipma'
plot(x, type = c("map", "path", "bottleneck"), scale = c("zoom", "full"), ...)
```

## Arguments

- x:

  A cipma object returned by
  [`cipma`](https://marekdejauj.github.io/cIPMA/reference/cipma.md).

- type:

  One of `"map"`, `"path"` (effect size across tolerances) or
  `"bottleneck"` (required level at the target level across tolerances).

- scale:

  One of `"zoom"` or `"full"`; used by the map.

- ...:

  Ignored.

## Value

A ggplot object.
