# Changelog

## cIPMA 0.2.0

- Quantile Necessary Condition Analysis (QNCA) of the ‘qnca’ package
  replaces Necessary Condition Analysis of the ‘NCA’ package as the
  necessity engine of
  [`cipma()`](https://marekdejauj.github.io/cIPMA/reference/cipma.md).
- Necessity is decided at the strict limit `pi = 1`, where the QNCA
  frontier equals the CE-FDH ceiling: the exact CE-FDH effect size and
  the QNCA permutation p-value at `pi = 1`, with the rule d \>= 0.10 and
  p \< 0.05.
- New argument `pi`: the tolerance of the bottleneck. The required level
  at the target level and the share of failing cases follow the chosen
  tolerance, and the effect size and p-value at the tolerance stand
  beside the strict-limit results.
- New tolerance path: the effect size and the required level at every
  tolerance of `pi_path`, by default 1.00 to 0.90 in steps of 0.01
  (`Tolerance_Path`, `plot(type = "path")`,
  `plot(type = "bottleneck")`).
- New
  [`cipma_diagnostics()`](https://marekdejauj.github.io/cIPMA/reference/cipma_diagnostics.md):
  the reference band, the deletion screen, the resolution table, the
  consistency probe and the addition screen of QNCA for every predictor.
- New [`summary()`](https://rdrr.io/r/base/summary.html) method with the
  detail levels `"brief"`, `"standard"` and `"full"`, and
  [`cipma_report()`](https://marekdejauj.github.io/cIPMA/reference/cipma_report.md)
  for a full Markdown report.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) gains
  `scale = "full"`, with performance on the whole 0-100 scale and
  importance from 0; `scale = "zoom"` keeps the axes of 0.1.6. All plots
  are in grey scale.
- The share of cases failing a bottleneck equals the percentile
  bottleneck table of the NCA package, the source of the failing shares
  of the original cIPMA (Hauff et al. 2024, Table 5): the scores are
  compared with the unrounded required level. Version 0.1.6 compared the
  scores with the required level rounded to three decimals, which added
  the case that sets the requirement whenever the rounding went up. On
  the technology acceptance data of Hauff et al. (2024) every required
  level at the target level 85 rounds down, and both versions reproduce
  Table 5 of the article. Required levels print with three decimals, as
  in the bottleneck table of NCA.
- Result columns: `NCA_d` becomes `d`, `NCA_p_value` becomes `p_value`,
  and `d_pi` and `p_pi` hold the results at the chosen tolerance. The
  bottleneck tables carry a `Status` column (`"required"`,
  `"not necessary"`, `"not attained"`), and `Bottleneck_full` is a long
  table across outcome levels.
- The argument `nca_rep` is deprecated; `B` sets the number of
  permutations, with default 1999.
- DESCRIPTION carries `Authors@R` and the issues address of the
  repository. The examples use the `mobi` data of ‘seminr’ and the
  technology acceptance data of Richter et al. (2023), read from
  Mendeley Data.

## cIPMA 0.1.6

- Combined importance-performance map analysis with PLS-SEM (‘seminr’)
  and Necessary Condition Analysis (‘NCA’).
