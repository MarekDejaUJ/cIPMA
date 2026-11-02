# cIPMA 0.2.0

- Quantile Necessary Condition Analysis (QNCA) of the 'qnca' package replaces Necessary Condition Analysis of the 'NCA' package as the necessity engine of `cipma()`.
- Necessity is decided at the strict limit `pi = 1`, where the QNCA frontier equals the CE-FDH ceiling: the exact CE-FDH effect size and the QNCA permutation p-value at `pi = 1`, with the rule d >= 0.10 and p < 0.05.
- New argument `pi`: the tolerance of the bottleneck. The required level at the target level and the share of failing cases follow the chosen tolerance, and the effect size and p-value at the tolerance stand beside the strict-limit results.
- New tolerance path: the effect size and the required level at every tolerance of `pi_path`, by default 1.00 to 0.90 in steps of 0.01 (`Tolerance_Path`, `plot(type = "path")`, `plot(type = "bottleneck")`).
- New `cipma_diagnostics()`: the reference band, the deletion screen, the resolution table, the consistency probe and the addition screen of QNCA for every predictor.
- New `summary()` method with the detail levels `"brief"`, `"standard"` and `"full"`, and `cipma_report()` for a full Markdown report.
- `plot()` gains `scale = "full"`, with performance on the whole 0-100 scale and importance from 0; `scale = "zoom"` keeps the axes of 0.1.6. All plots are in grey scale.
- The share of cases failing a bottleneck uses the unrounded required level. Version 0.1.6 compared the scores with the required level rounded to three decimals, which counted the case that sets the requirement as failing whenever the rounding went up.
- Result columns: `NCA_d` becomes `d`, `NCA_p_value` becomes `p_value`, and `d_pi` and `p_pi` hold the results at the chosen tolerance. The bottleneck tables carry a `Status` column (`"required"`, `"not necessary"`, `"not attained"`), and `Bottleneck_full` is a long table across outcome levels.
- The argument `nca_rep` is deprecated; `B` sets the number of permutations, with default 1999.
- DESCRIPTION carries `Authors@R` and the issues address of the repository. The examples use the `mobi` data of 'seminr' and the technology acceptance data of Richter et al. (2023), read from Mendeley Data.

# cIPMA 0.1.6

- Combined importance-performance map analysis with PLS-SEM ('seminr') and Necessary Condition Analysis ('NCA').
