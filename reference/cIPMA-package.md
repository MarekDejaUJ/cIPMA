# cIPMA: Combined Importance-Performance Map Analysis

The cIPMA package combines the importance-performance map analysis of a
PLS-SEM model, estimated with seminr, with Quantile Necessary Condition
Analysis of each predictor, computed with qnca.
[`cipma`](https://marekdejauj.github.io/cIPMA/reference/cipma.md) runs
the analysis; [`summary()`](https://rdrr.io/r/base/summary.html) reports
the results at three levels of detail;
[`cipma_diagnostics`](https://marekdejauj.github.io/cIPMA/reference/cipma_diagnostics.md)
runs the QNCA diagnostics;
[`cipma_report`](https://marekdejauj.github.io/cIPMA/reference/cipma_report.md)
writes a full report;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
combined map and the tolerance paths.

## References

Hauff, S., Richter, N. F., Sarstedt, M., & Ringle, C. M. (2024).
Importance and performance in PLS-SEM and NCA: Introducing the combined
importance-performance map analysis (cIPMA). Journal of Retailing and
Consumer Services, 78, 103723.
[doi:10.1016/j.jretconser.2024.103723](https://doi.org/10.1016/j.jretconser.2024.103723)

Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and
methodology of "necessary but not sufficient" causality. Organizational
Research Methods, 19(1), 10-52.
[doi:10.1177/1094428115584005](https://doi.org/10.1177/1094428115584005)
