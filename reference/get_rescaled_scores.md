# Get Rescaled Latent Variable Scores

Extracts latent variable scores from a seminr model and rescales them to
0-100.

## Usage

``` r
get_rescaled_scores(model, scales)
```

## Arguments

- model:

  A seminr model object.

- scales:

  A named list where keys are indicator names and values are c(min,
  max).

## Value

A data frame of rescaled latent variable scores.
