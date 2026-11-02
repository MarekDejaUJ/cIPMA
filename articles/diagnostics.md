# QNCA diagnostics in cIPMA

[`cipma_diagnostics()`](https://marekdejauj.github.io/cIPMA/reference/cipma_diagnostics.md)
runs the diagnostics of Quantile Necessary Condition Analysis (QNCA;
Deja 2026) on every predictor of a combined importance-performance map
analysis (cIPMA), with the rescaled latent variable scores and the
empirical scope of each pair of predictor and target construct. The
scope is the rectangle between the observed minimum and maximum of both
scores. The tolerance `pi` in (0, 1\] lets the share `1 - pi` of the
successful cases, the cases that attain a target level, sit below the
requirement of the predictor, and the effect size `d_pi` is the share of
the scope in the empty zone to the left of the QNCA frontier at `pi`. At
`pi = 1` the QNCA frontier equals the CE-FDH ceiling of Necessary
Condition Analysis (NCA; Dul 2016).

The reference band compares the observed effect size with a benchmark
that keeps the dependence of the pair, the deletion screen and the
addition screen measure the influence of single cases, the resolution
table gives the rank behind each bottleneck, and the consistency probe
follows the effect size across subsample sizes.

The example uses the `mobi` data of `seminr`, customers of a mobile
phone provider in the European customer satisfaction index model, with
loyalty at level 80 of 100 as the target and the tolerance `pi = 0.95`
for the bottleneck.

``` r

library(cIPMA)
library(seminr)
mm <- constructs(
  composite("Image", multi_items("IMAG", 1:5)),
  composite("Expectation", multi_items("CUEX", 1:3)),
  composite("Quality", multi_items("PERQ", 1:7)),
  composite("Value", multi_items("PERV", 1:2)),
  composite("Satisfaction", multi_items("CUSA", 1:3)),
  composite("Complaints", single_item("CUSCO")),
  composite("Loyalty", multi_items("CUSL", 1:3))
)
sm <- relationships(
  paths(from = "Image", to = c("Expectation", "Satisfaction", "Loyalty")),
  paths(from = "Expectation", to = c("Quality", "Value", "Satisfaction")),
  paths(from = "Quality", to = c("Value", "Satisfaction")),
  paths(from = "Value", to = "Satisfaction"),
  paths(from = "Satisfaction", to = c("Complaints", "Loyalty")),
  paths(from = "Complaints", to = "Loyalty")
)
pls <- estimate_pls(mobi, mm, sm)
boot <- bootstrap_model(pls, nboot = 200, cores = 1, seed = 123)
scales <- setNames(rep(list(c(1, 10)), ncol(mobi)), names(mobi))
res <- cipma(boot, "Loyalty", mobi, scales, target_level = 80, pi = 0.95,
             B = 499, seed = 1, pls_model = pls)
res
#> cIPMA: target Loyalty at level 80 (n = 250, tolerance pi = 0.95)
#> 
#>     Construct Importance Performance     d      p Necessary Required Failing
#>         Image      0.399        73.5 0.151 0.0140       yes   63.566   21.6%
#>   Expectation      0.218        72.5 0.122 0.0100       yes   56.419   10.0%
#>       Quality      0.323        75.0 0.151 0.0040       yes   62.839   16.4%
#>         Value      0.101        62.3 0.016 0.5600        no   44.444   14.8%
#>  Satisfaction      0.521        72.5 0.140 0.0020       yes   67.209   29.2%
#>    Complaints      0.067        67.4 0.016 0.3720        no   44.444   10.8%
#> 
#> Necessary: CE-FDH d >= 0.10 and p < 0.05 at pi = 1.00. Required and Failing:
#> bottleneck at pi = 0.95 (NN: no requirement; NA: level not attained).
#> summary(x) adds the readings; summary(x, detail = "full") adds the tolerance path and
#> the bottleneck table; cipma_diagnostics(x) runs the QNCA diagnostics.
```

The example draws 200 bootstrap resamples and 499 permutations, which
keeps the computation short. The number of bootstrap resamples is a
choice of the analyst in
[`bootstrap_model()`](https://rdrr.io/pkg/seminr/man/bootstrap_model.html),
and [`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md)
draws 1999 permutations by default.

## Running the diagnostics

The argument `points` of
[`cipma_diagnostics()`](https://marekdejauj.github.io/cIPMA/reference/cipma_diagnostics.md)
declares the position of the cases that the addition screen adds. The
position `c(30, 95)` describes an efficient success, with a predictor
score of 30 and loyalty of 95, and the position lies inside the
empirical scope of every pair.

``` r

round(sapply(res$Data_Rescaled, range), 1)
#>      Image Expectation Quality Value Satisfaction Complaints Loyalty
#> [1,]  26.9        26.7    24.4     0         27.6          0     0.9
#> [2,] 100.0       100.0   100.0   100        100.0        100   100.0
```

The highest of the lowest predictor scores is 27.6, for satisfaction,
and loyalty reaches 100, so the added cases at predictor score 30 and
loyalty 95 fall inside every scope.

``` r

dg <- cipma_diagnostics(res, points = c(30, 95), seed = 1)
dg
#> QNCA diagnostics of the cIPMA predictors of Loyalty (tolerances 1.00, 0.95, 0.90)
#> 
#> Reference band (199 benchmark datasets): observed d, upper band limit, excess
#>     Construct   pi d_obs upper excess p_null Above_band
#>         Image 1.00 0.152 0.271 -0.118  0.780         no
#>         Image 0.95 0.404 0.466 -0.062  0.730         no
#>         Image 0.90 0.496 0.517 -0.021  0.460         no
#>   Expectation 1.00 0.123 0.225 -0.102  0.630         no
#>   Expectation 0.95 0.363 0.419 -0.056  0.850         no
#>   Expectation 0.90 0.446 0.473 -0.027  0.495         no
#>       Quality 1.00 0.152 0.286 -0.134  0.840         no
#>       Quality 0.95 0.420 0.487 -0.066  0.720         no
#>       Quality 0.90 0.526 0.557 -0.030  0.405         no
#>         Value 1.00 0.014 0.200 -0.187  1.000         no
#>         Value 0.95 0.358 0.418 -0.060  0.830         no
#>         Value 0.90 0.443 0.489 -0.047  0.640         no
#>  Satisfaction 1.00 0.141 0.258 -0.117  0.780         no
#>  Satisfaction 0.95 0.402 0.428 -0.026  0.315         no
#>  Satisfaction 0.90 0.471 0.494 -0.023  0.290         no
#>    Complaints 1.00 0.014 0.109 -0.095  0.965         no
#>    Complaints 0.95 0.316 0.406 -0.090  0.405         no
#>    Complaints 0.90 0.440 0.497 -0.057  0.575         no
#> 
#> Deletion screen (k = 1): most influential deletion and its relative change of d per tolerance
#>     Construct Cases       Relative_change Flagged_rows
#>         Image    28 0.191 / 0.044 / 0.007           17
#>   Expectation   117 0.302 / 0.010 / 0.006            5
#>       Quality   133 0.258 / 0.057 / 0.012           18
#>         Value   242 2.403 / 0.016 / 0.028           56
#>  Satisfaction   188 0.194 / 0.009 / 0.008            6
#>    Complaints   109 0.000 / 0.036 / 0.010           26
#> 
#> Resolution at pi = 0.95 for outcome targets up to 80
#>     Construct Targets Inherited Smallest_resistance Smallest_attaining_set
#>         Image      40         6                   5                    104
#>   Expectation      40        11                   5                    104
#>       Quality      40        10                   5                    104
#>         Value      40         8                   5                    104
#>  Satisfaction      40         2                   5                    104
#>    Complaints      40         3                   5                    104
#> 
#> Consistency probe: divergence of the drift slopes and the slopes at the two tolerances
#>     Construct Divergence              Drift
#>         Image     10.008 10.26297 / 0.25493
#>   Expectation      8.591  9.20613 / 0.61513
#>       Quality      6.953 7.89470 / -0.94169
#>         Value     10.615 11.37308 / 0.75806
#>  Satisfaction     10.156 10.55362 / 0.39722
#>    Complaints      6.167 6.80132 / -0.63434
#> 
#> An empty zone above the band exceeds the dependence-matched benchmark; the band is a
#> proposed diagnostic with operating characteristics still to be established.
```

The diagnostics use the tolerances 1.00, 0.95 and 0.90 and a grid of 50
outcome levels, 199 benchmark datasets for the band and 50 subsamples
per size for the probe, the defaults of
[`cipma_diagnostics()`](https://marekdejauj.github.io/cIPMA/reference/cipma_diagnostics.md).
On the grid, the effect size at `pi = 1` differs from the exact CE-FDH
effect size of
[`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md) by
the error of the grid: 0.152 against 0.151 for image, and 0.014 against
0.016 for value.

## Reference band

The reference band draws benchmark datasets from a Gaussian copula that
keeps the empirical marginals of the predictor and the target construct
and the normal-score rank correlation of the pair, so the scatter of
each benchmark dataset follows a smooth monotone dependence, free of an
attainability floor (Deja 2026). The full QNCA estimator runs on every
benchmark dataset at every tolerance on the same grid and scope. The
column `upper` gives the 97.5th percentile of the benchmark effect
sizes, `excess` the observed effect size minus the upper limit, and
`p_null` the share of benchmark effect sizes at least as large as the
observed effect size, with one added to the count and to the number of
datasets. `Above_band` reads yes when the excess is positive.

Every excess in the example is negative, from -0.187 for value at
`pi = 1.00` to -0.021 for image at `pi = 0.90`, and every predictor
stays within the band at the three tolerances. Within the band, the
empty zone keeps the size that the marginals and the rank dependence of
the pair reproduce alone. Image, expectation, quality and satisfaction
meet the necessity rule of
[`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md) with
the permutation test, which compares the effect size with independent
shuffles of the target construct, and the band withholds support for
necessity beyond the dependence-matched benchmark and leaves the
necessity hypotheses open.

The band targets spuriousness, an empty zone produced by dependence and
skewness, a failure mode distinct from outlier-fragility. The band
benchmarks the effect size and leaves the QNCA estimator, the effect
size and the bottleneck as defined. The reference band is a proposed
diagnostic with operating characteristics still to be established (Deja
2026). The Gaussian copula reproduces a smooth monotone dependence and
leaves arbitrary tail dependence, latent classes and several sufficient
pathways outside the benchmark, and the band holds pointwise at each
tolerance.

## Deletion screen

The deletion screen removes each case in turn (`k = 1`) and recomputes
the effect size on the scope of the full sample at each tolerance. The
printed table gives, for each predictor, the case whose deletion changes
the effect size most, the relative change of `d_pi` at 1.00, 0.95 and
0.90 after the deletion of the case, and the number of deletions with a
relative change of at least 0.01 at some tolerance (`Flagged_rows`). The
element `outliers` of each predictor holds the full screen.

``` r

head(dg$results$Image$outliers$table, 3)
#>   cases flagged d_without_1  dif_rel_1 d_without_0.95 dif_rel_0.95
#> 1    28    TRUE   0.1814358 0.19073843      0.4219372   0.04423097
#> 2   242    TRUE   0.1595864 0.04734373      0.4219372   0.04423097
#> 3   133    TRUE   0.1596182 0.04755270      0.4219265   0.04420449
#>   d_without_0.9 dif_rel_0.9
#> 1     0.4997136 0.007130366
#> 2     0.5005153 0.008746165
#> 3     0.5021290 0.011998429
```

For image, the deletion of case 28 raises the effect size at `pi = 1`
from 0.152 to 0.181, a relative change of 0.191, with relative changes
of 0.044 at `pi = 0.95` and 0.007 at `pi = 0.90`. The strict boundary of
image rests on one customer, whereas the tolerant members change little
after the same deletion, and the deletion screen reports the
outlier-fragility of the strict member that the tolerance path, the
effect size from `pi = 1.00` to `pi = 0.90`, also shows (Deja 2026).
Value shows the largest relative change, 2.403 for case 242 at
`pi = 1.00`, because the strict effect size of value lies close to zero.

## Resolution table

The resolution table follows the bottleneck at the tolerance of the
cIPMA, `pi = 0.95`, target by target on the grid of 50 outcome levels
(Deja 2026). For each target the table gives the number `k` of
successful cases, the selected rank `max(1, ceiling(k * (1 - pi)))`, the
raw order statistic, the requirement after the running maximum across
targets, the flag of a requirement inherited from a lower target, the
number and share of successful cases below the requirement, and the
resistance, the rank minus 1. The printed table counts, for the targets
up to the target level, the targets, the inherited requirements, the
smallest resistance and the smallest attaining set.

``` r

r <- dg$results$Image$resolution
r <- r[r$y >= 65 & r$y <= 80, ]
data.frame(y = round(r$y, 2), k = r$k, rank = r$rank, raw = round(r$raw, 3),
           fitted = round(r$fitted, 3), inherited = r$inherited, below = r$below,
           share = round(r$exception_share, 3), resistance = r$resistance)
#>       y   k rank    raw fitted inherited below share resistance
#> 1 65.61 162    9 59.357 59.357     FALSE     8 0.049          8
#> 2 67.63 153    8 57.996 59.357      TRUE     8 0.052          7
#> 3 69.66 143    8 57.996 59.357      TRUE     8 0.056          7
#> 4 71.68 137    7 57.996 59.357      TRUE     7 0.051          6
#> 5 73.70 134    7 59.466 59.466     FALSE     6 0.045          6
#> 6 75.73 129    7 59.627 59.627     FALSE     6 0.047          6
#> 7 77.75 112    6 61.793 61.793     FALSE     5 0.045          5
#> 8 79.77 104    6 63.566 63.566     FALSE     5 0.048          5
```

At the highest target up to 80, loyalty of 79.77, 104 customers attain
the target, the rank is 6 and the requirement of image is 63.566 with
resistance 5: the raw requirement absorbs five added efficient cases and
stays at or above the lowest image score of the original successful
customers (Deja 2026). At loyalty of 67.63, 69.66 and 71.68 the
requirement of 59.357 is inherited from a lower target, and the share of
successful customers below the requirement, 0.052, 0.056 and 0.051,
exceeds `1 - pi`: an inherited requirement can exceed the tolerance when
the exceptions admitted at a lower target concentrate among the cases
that attain the higher target. Over the 40 targets up to 80, image
inherits the requirement at 6 targets, and every predictor keeps a
resistance of at least 5 and an attaining set of at least 104 customers.

## Consistency probe

The consistency probe draws 50 subsamples without replacement at each of
four sizes, a quarter, a half, three quarters and the whole of the
sample, and averages `d_pi` at `pi = 1` and at the lowest tolerance of
the grid, `pi = 0.90`, on the scope of the full sample. The drift slope
is the least-squares slope of the mean effect size on the inverse
subsample size, and the divergence is the absolute drift slope at
`pi = 1` minus the absolute drift slope at `pi = 0.90`.

``` r

dg$results$Image$probe$summary
#>   size  pi    mean_d       sd_d
#> 1   62 1.0 0.2800864 0.07745286
#> 2  125 1.0 0.2170544 0.04236041
#> 3  186 1.0 0.1718910 0.01882249
#> 4  250 1.0 0.1523725 0.00000000
#> 5   62 0.9 0.4979799 0.03652387
#> 6  125 0.9 0.5000494 0.01822201
#> 7  186 0.9 0.4928405 0.01029921
#> 8  250 0.9 0.4961757 0.00000000
```

The mean strict effect size of image falls from 0.280 in subsamples of
62 customers to 0.152 in the whole sample, whereas the mean effect size
at `pi = 0.90` stays between 0.493 and 0.500. The strict effect size,
built on the minimum of each attaining set, depends on the sample size,
and the tolerant effect size stays stable; the drift slopes of image are
10.263 and 0.255, with a divergence of 10.008. Every predictor shows a
positive divergence, from 6.167 for complaints to 10.615 for value. The
probe is a descriptive check of subsample stability, and hard bounds of
the support, for example the limits of a rescaled 0–100 scale, and
skewed marginals can confound the probe.

## Addition screen

With `points`, the addition screen adds 0 to 5 cases at the declared
position and recomputes `d_pi` on the scope of the original sample (Deja
2026). The addition screen stays outside the printed tables, and the
element `sensitivity` of each predictor holds the effect sizes.

``` r

s <- dg$results$Image$sensitivity
d_added <- round(s$d, 3)
dimnames(d_added) <- list(added = s$counts, pi = format(s$pi))
d_added
#>      pi
#> added  1.00  0.95  0.90
#>     0 0.152 0.404 0.496
#>     1 0.049 0.392 0.486
#>     2 0.049 0.366 0.483
#>     3 0.049 0.329 0.477
#>     4 0.049 0.311 0.468
#>     5 0.049 0.303 0.455
```

One added case at image 30 and loyalty 95 lowers the strict effect size
of image from 0.152 to 0.049, because the added case attains every
target up to 95 and lowers the strict requirement of each of the targets
to at most 30. At `pi = 0.95` one added case lowers `d_pi` from 0.404 to
0.392 and five added cases to 0.303, and at `pi = 0.90` five added cases
lower `d_pi` from 0.496 to 0.455. At a target with rank `h` above the
number `c` of added cases, the requirement stays at or above the order
statistic of rank `h - c` among the original successful customers. Each
added case therefore lowers the requirement by at most one spacing
between the image scores of the original successful customers, and the
loss of the tolerant members grows step by step with the number of added
cases (Deja 2026).

## Report

[`cipma_report()`](https://marekdejauj.github.io/cIPMA/reference/cipma_report.md)
writes the full summary of the cIPMA as a Markdown report, and adds the
tables of the diagnostics with the argument `diagnostics`.

``` r

report_file <- tempfile(fileext = ".md")
report <- cipma_report(res, file = report_file, diagnostics = dg)
length(report)
#> [1] 144
grep("^#", report, value = TRUE)
#> [1] "# cIPMA report: Loyalty at level 80"                       
#> [2] "## Combined map"                                           
#> [3] "## PLS-SEM significance and effect at the chosen tolerance"
#> [4] "## Readings"                                               
#> [5] "## Tolerance path: effect size d"                          
#> [6] "## Tolerance path: required level for Loyalty >= 80"       
#> [7] "## Bottleneck table at pi = 0.95"                          
#> [8] "## QNCA diagnostics"
```

The report of the example holds 144 lines, with the combined map, the
PLS-SEM significance, the readings, the two tolerance paths, the
bottleneck table and the QNCA diagnostics.

## References

- Deja, M. (2026). Quantile Necessary Condition Analysis: A
  Tolerance-Parameterized Generalization of Necessary Condition
  Analysis. *Quality & Quantity* (accepted for publication).
- Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and
  methodology of “necessary but not sufficient” causality.
  *Organizational Research Methods*, 19(1), 10–52.
  <https://doi.org/10.1177/1094428115584005>
