#' Print cIPMA Results
#'
#' @param x A cIPMA object.
#' @param ... Ignored.
#'
#' @export
print.cipma <- function(x, ...) {
  cat("=== Combined Importance-Performance Map Analysis (cIPMA) ===\n")
  cat("Target Construct:", x$Dependent, "\n")
  cat("Desired Outcome Level:", x$Target_Level, "\n\n")

  cat("--- Table A: PLS-SEM & NCA Statistics ---\n")
  print(x$cIPMA_results[, c("Construct", "Importance", "Performance",
                            "PLS_p_value", "NCA_d", "NCA_p_value", "Is_Necessary")])

  cat("\n--- Table B: Bottleneck Summary ---\n")
  print(x$Bottleneck_summary)
}


utils::globalVariables(c("Importance", "Performance", "PlotSize", "ColorGroup", "Label", "Fail_Percentage", "Is_Necessary", "Construct"))
