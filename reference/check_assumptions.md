# Check cIPMA Assumptions and Model Diagnostics

Check cIPMA Assumptions and Model Diagnostics

## Usage

``` r
check_assumptions(x, vif_threshold = 5, htmt_threshold = 0.9)
```

## Arguments

- x:

  A cipma object or a seminr model object.

- vif_threshold:

  Numeric. VIF threshold for multicollinearity. Default is 5.

- htmt_threshold:

  Numeric. HTMT threshold for discriminant validity. Default is 0.90.

## Value

A list containing assumption checks and diagnostics.
