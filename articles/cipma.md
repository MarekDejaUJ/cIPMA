# The cIPMA procedure

The combined importance-performance map analysis (cIPMA) of Hauff,
Richter, Sarstedt and Ringle (2024) places each predictor of a target
construct on a map of importance against performance and marks the
predictors that are necessary for a target level of the construct.
Importance is the total effect of the predictor on the target construct
in a model of partial least squares structural equation modelling
(PLS-SEM), and performance is the mean of the latent variable scores of
the predictor, rescaled to 0–100. Necessity comes from Necessary
Condition Analysis (NCA), which finds the level of a predictor required
for a target level of the outcome (Dul 2016).

The cIPMA package computes necessity with Quantile Necessary Condition
Analysis (QNCA; Deja 2026). At the tolerance `pi = 1` the QNCA frontier
equals the CE-FDH ceiling of NCA, the ceiling of the published cIPMA,
and the package decides necessity at `pi = 1`. A tolerance `pi < 1` lets
the share `1 - pi` of the cases that attain the target level sit below
the reported requirement of a predictor.

## Steps

Hauff et al. (2024, Fig. 3) organise the cIPMA in eight steps. The
analyst carries out the first three steps with the theory of the study
and the `seminr` package, and the cIPMA package supports the last five:

1.  Specify the research objective and the theoretical background, with
    hypotheses in the sufficiency logic of PLS-SEM and in the necessity
    logic of NCA.
2.  Prepare and check the data, with the indicators of each construct on
    one common scale, for example 1 to 5.
3.  Run the PLS-SEM analysis with
    [`estimate_pls()`](https://rdrr.io/pkg/seminr/man/estimate_pls.html)
    and
    [`bootstrap_model()`](https://rdrr.io/pkg/seminr/man/bootstrap_model.html)
    of `seminr`.
4.  Evaluate the reliability and validity of the measurement models and
    check that every outer weight is positive, with
    [`check_assumptions()`](https://marekdejauj.github.io/cIPMA/reference/check_assumptions.md).
5.  Rescale the latent variable scores to 0–100 with the theoretical
    minimum and maximum of each indicator, inside
    [`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md).
6.  Run the necessity analysis on the rescaled scores with the empirical
    scope and read the bottleneck at the target level, inside
    [`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md).
7.  Evaluate the PLS-SEM significance of the total effects and the
    effect size and p-value of necessity, with
    [`print()`](https://rdrr.io/r/base/print.html) and
    [`summary()`](https://rdrr.io/r/base/summary.html).
8.  Draw the combined map with
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Data and model

The example uses the `mobi` data of `seminr`: customers of a mobile
phone provider, with indicators on a scale from 1 to 10, in the European
customer satisfaction index model.

``` r

library(cIPMA)
library(seminr)
dim(mobi)
#> [1] 250  24
range(unlist(mobi))
#> [1]  1 10
```

The 250 customers answer 24 indicators. Loyalty is the target construct,
and image, expectation, quality, value, satisfaction and complaints act
on loyalty directly or through other constructs.

``` r

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
```

The example draws 200 bootstrap resamples and, inside
[`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md), 499
permutations, which keeps the computation short. The number of bootstrap
resamples is a choice of the analyst in
[`bootstrap_model()`](https://rdrr.io/pkg/seminr/man/bootstrap_model.html),
and [`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md)
draws 1999 permutations by default.

## Measurement models

Step 4 checks the measurement models with the criteria of PLS-SEM and
adds the requirement of the importance-performance map analysis (IPMA):
every outer weight is positive, which keeps the rescaled scores inside
0–100 (Hauff et al. 2024).

``` r

check_assumptions(pls)
#> === cIPMA Assumption Checks ===
#> 
#> 1. POSITIVE OUTER WEIGHTS
#>    [PASS] All outer weights are positive.
#> 
#> 2. RELIABILITY (alpha >= 0.70, rhoC >= 0.70, AVE >= 0.50)
#>              alpha  rhoA  rhoC   AVE
#> Image        0.723 0.740 0.819 0.478
#> Expectation  0.452 0.462 0.733 0.480
#> Quality      0.877 0.884 0.905 0.577
#> Value        0.824 0.855 0.918 0.848
#> Satisfaction 0.779 0.789 0.871 0.693
#> Complaints   1.000 1.000 1.000 1.000
#> Loyalty      0.472 0.746 0.722 0.517
#> 
#> Alpha, rhoA, and rhoC should exceed 0.7 while AVE should exceed 0.5
#>    [WARN] Check constructs: Image, Expectation, Loyalty 
#> 
#> 3. DISCRIMINANT VALIDITY (HTMT < 0.9 )
#>              Image Expectation Quality Value Satisfaction Complaints Loyalty
#> Image            .           .       .     .            .          .       .
#> Expectation  0.888           .       .     .            .          .       .
#> Quality      0.929       0.878       .     .            .          .       .
#> Value        0.652       0.589   0.673     .            .          .       .
#> Satisfaction 0.910       0.865   0.954 0.741            .          .       .
#> Complaints   0.545       0.383   0.564 0.387        0.588          .       .
#> Loyalty      0.867       0.770   0.759 0.797        0.957      0.561       .
#>    [FAIL] HTMT exceeded for:
#>           Quality <-> Image 
#>           Satisfaction <-> Image 
#>           Satisfaction <-> Quality 
#>           Loyalty <-> Satisfaction 
#> 
#> 4. MULTICOLLINEARITY (VIF < 5 )
#>    Expectation :
#>     Image = NA 
#>    Satisfaction :
#>     Image = 2.37, Expectation = 1.49, Quality = 2.84, Value = 1.55 
#>    Loyalty :
#>     Image = 1.99, Satisfaction = 2.13, Complaints = 1.43 
#>    Quality :
#>     Expectation = NA 
#>    Value :
#>     Expectation = 1.45, Quality = 1.45 
#>    Complaints :
#>     Satisfaction = NA 
#>    [PASS] No multicollinearity issues detected.
#> 
#> 5. SAMPLE SIZE
#>    N = 250 (recommended >= 70 )
#>    [PASS] Sample size adequate.
```

Every outer weight of the model is positive, and the variance inflation
factors of the structural model stay below 5. The reliability check
flags image, expectation and loyalty, and four heterotrait-monotrait
ratios exceed 0.90. Hauff et al. (2024) ask for improvements of a
measurement model that fails the assessment criteria, and the example
keeps the model unchanged as a demonstration of the procedure.

## Rescaled scores, importance and performance

[`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md)
computes the rescaled latent variable scores of Hauff et al. (2024). The
standardised outer weights, divided by the standard deviation of each
indicator, are normalised to sum to one within each construct, and the
weighted sum of the indicators, rescaled with the theoretical minimum
and maximum of the indicators given in the named list `scales`, gives
the score of each case on 0–100. Importance is the total effect of each
predictor on the target construct, the direct effect plus every indirect
effect through mediating constructs, and performance is the mean
rescaled score.

The target level of loyalty is 80 on the 0–100 scale. Hauff et
al. (2024) derive the target level from expert knowledge, from the
standards of the field or from the aims of the project. The bottleneck
uses the tolerance `pi = 0.95`.

``` r

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

Satisfaction has the largest importance, 0.521, followed by image with
0.399, and complaints has the smallest, 0.067. Performance ranges from
62.3 for value to 75.0 for quality.

## Necessity and the bottleneck

For each predictor,
[`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md)
computes the CE-FDH effect size `d` of NCA on the empirical scope of the
rescaled pair: the area of the empty zone in the upper-left corner of
the scatter of predictor and target construct, divided by the area of
the scope, the rectangle between the observed minimum and maximum of
both scores. The permutation test shuffles the target construct against
the predictor `B` times, recomputes the effect size at `pi = 1` on a
grid of 100 outcome levels after each shuffle, and reports the share of
permuted effect sizes at least as large as the observed effect size,
with one added to the count and to the number of permutations. A
predictor is necessary when `d` reaches 0.10 and the p-value stays below
0.05, the rule of Hauff et al. (2024) with the effect-size landmark of
Dul (2016).

Image, expectation, quality and satisfaction meet the rule, with effect
sizes from 0.122 to 0.151 and p-values from 0.0020 to 0.0140. A p-value
of 0.0020 is the smallest value that 499 permutations report. Value and
complaints miss the rule, each with `d` = 0.016, and with p-values of
0.5600 and 0.3720.

The column `Required` gives the bottleneck: the level of the predictor
required for loyalty of at least 80. At `pi = 0.95` the requirement is
the 5% quantile of the predictor among the successful customers, the
customers with loyalty of at least 80, so the most efficient 5% of the
successful customers may sit below the requirement. The column `Failing`
gives the share of all customers below the required level, which Hauff
et al. (2024) read as the performance of the predictor in the necessity
logic. Satisfaction carries the highest requirement, 67.209, and the
largest share of customers below the requirement, 29.2%. Value and
complaints also carry a requirement of 44.444 at `pi = 0.95`, and the
decision at `pi = 1` keeps both outside the necessary predictors.

## Results at three levels of detail

`summary(res)` adds the bootstrap p-values of the total effects, the
effect size and p-value at the chosen tolerance, the number of customers
below each requirement and one reading per predictor.

``` r

summary(res)
#> Combined importance-performance map analysis
#> Target construct: Loyalty, target level: 80, cases: 250
#> Tolerance of the bottleneck: pi = 0.95; permutations: 499; outcome levels: 400 (test: 100); qnca 0.4.0
#> 
#> Combined map
#>     Construct Importance Performance     d      p Necessary Required Failing
#>         Image      0.399        73.5 0.151 0.0140       yes   63.566   21.6%
#>   Expectation      0.218        72.5 0.122 0.0100       yes   56.419   10.0%
#>       Quality      0.323        75.0 0.151 0.0040       yes   62.839   16.4%
#>         Value      0.101        62.3 0.016 0.5600        no   44.444   14.8%
#>  Satisfaction      0.521        72.5 0.140 0.0020       yes   67.209   29.2%
#>    Complaints      0.067        67.4 0.016 0.3720        no   44.444   10.8%
#> 
#> PLS-SEM significance, effect at the chosen tolerance, failing cases
#>     Construct   PLS_p d_0.95 p_0.95 Failing_cases
#>         Image <0.0001  0.404 0.0020            54
#>   Expectation <0.0001  0.363 0.0580            25
#>       Quality <0.0001  0.421 0.0020            41
#>         Value  0.0143  0.359 0.0020            37
#>  Satisfaction <0.0001  0.405 0.0020            73
#>    Complaints  0.2906  0.316 0.0020            27
#> 
#> Readings
#> - Image is necessary at the strict limit (d = 0.151, p = 0.0140); the required level of Image for Loyalty >= 80 at pi = 0.95 is 63.566, with 21.6% of the cases below the required level.
#> - Expectation is necessary at the strict limit (d = 0.122, p = 0.0100); the required level of Expectation for Loyalty >= 80 at pi = 0.95 is 56.419, with 10.0% of the cases below the required level.
#> - Quality is necessary at the strict limit (d = 0.151, p = 0.0040); the required level of Quality for Loyalty >= 80 at pi = 0.95 is 62.839, with 16.4% of the cases below the required level.
#> - Value is not necessary at the strict limit (d = 0.016, p = 0.5600); the required level of Value for Loyalty >= 80 at pi = 0.95 is 44.444, with 14.8% of the cases below the required level.
#> - Satisfaction is necessary at the strict limit (d = 0.140, p = 0.0020); the required level of Satisfaction for Loyalty >= 80 at pi = 0.95 is 67.209, with 29.2% of the cases below the required level.
#> - Complaints is not necessary at the strict limit (d = 0.016, p = 0.3720); the required level of Complaints for Loyalty >= 80 at pi = 0.95 is 44.444, with 10.8% of the cases below the required level.
#> 
#> Effect sizes at pi = 0.95 are relative to the tolerance and rise as pi falls;
#> the NCA landmark of 0.10 applies to the CE-FDH effect size at pi = 1.00.
```

The total effects of image, expectation, quality and satisfaction have
bootstrap p-values below 0.0001, the total effect of value 0.0143 and
the total effect of complaints 0.2906. The effect sizes at `pi = 0.95`
range from 0.316 to 0.421 and exceed the CE-FDH effect sizes of every
predictor, because the QNCA frontier at a tolerance below 1 lies on or
above the CE-FDH ceiling. Effect sizes at `pi < 1` are relative to the
tolerance, and the NCA landmark of 0.10 applies to the CE-FDH effect
size at `pi = 1`.

`summary(res, detail = "brief")` gives the combined table with the
settings, and `summary(res, detail = "full")` adds the tolerance path of
the effect size and of the required level from `pi = 1.00` to
`pi = 0.90` and the bottleneck table across outcome levels.

## The combined map

`plot(res)` draws the combined map of step 8: importance on the
horizontal axis, performance on the vertical axis, and dashed lines at
the mean importance and the mean performance of the predictors. White
bubbles mark the predictors necessary at `pi = 1`, with a size given by
the share of customers below the required level, and black dots mark the
other predictors.

``` r

plot(res)
```

![Combined map of the predictors of loyalty at level 80 for the mobi
data. Importance is the total effect on loyalty, and performance is the
mean rescaled score on 0–100. White bubbles mark the predictors
necessary at pi = 1, with a size given by the share of customers below
the required level at pi = 0.95; black dots mark value and
complaints.](cipma_files/figure-html/procedure-map-1.png)

Combined map of the predictors of loyalty at level 80 for the mobi data.
Importance is the total effect on loyalty, and performance is the mean
rescaled score on 0–100. White bubbles mark the predictors necessary at
pi = 1, with a size given by the share of customers below the required
level at pi = 0.95; black dots mark value and complaints.

Satisfaction, image and quality lie in the quadrant of high importance
and high performance, and expectation lies in the quadrant of low
importance and high performance. Value and complaints lie in the
quadrant of low importance and low performance.

## References

- Deja, M. (2026). Quantile Necessary Condition Analysis: A
  Tolerance-Parameterized Generalization of Necessary Condition
  Analysis. *Quality & Quantity* (accepted for publication).
- Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and
  methodology of “necessary but not sufficient” causality.
  *Organizational Research Methods*, 19(1), 10–52.
  <https://doi.org/10.1177/1094428115584005>
- Hauff, S., Richter, N. F., Sarstedt, M., & Ringle, C. M. (2024).
  Importance and performance in PLS-SEM and NCA: Introducing the
  combined importance-performance map analysis (cIPMA). *Journal of
  Retailing and Consumer Services*, 78, 103723.
  <https://doi.org/10.1016/j.jretconser.2024.103723>
