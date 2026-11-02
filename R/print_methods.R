# Formatting helpers shared by print, summary and report.

fmt_num <- function(x, digits = 3) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- "NA"
  out
}

fmt_p <- function(p) {
  out <- formatC(p, format = "f", digits = 4)
  out[!is.na(p) & p < 1e-4] <- "<0.0001"
  out[is.na(p)] <- "NA"
  out
}

fmt_yes <- function(x) {
  out <- ifelse(x, "yes", "no")
  out[is.na(x)] <- "NA"
  out
}

fmt_required <- function(required, status, digits = 2) {
  out <- fmt_num(required, digits)
  out[status == "not necessary"] <- "NN"
  out[status == "not attained"] <- "NA"
  out
}

# The combined table of a cipma object, formatted for reading.
cipma_main_table <- function(x, digits = 3) {
  r <- x$cIPMA_results
  b <- x$Bottleneck_summary[match(r$Construct, x$Bottleneck_summary$Construct), ]
  data.frame(
    Construct = r$Construct,
    Importance = fmt_num(r$Importance, digits),
    Performance = fmt_num(r$Performance, 1),
    d = fmt_num(r$d, digits),
    p = fmt_p(r$p_value),
    Necessary = fmt_yes(r$Is_Necessary),
    Required = fmt_required(b$Required_Level, b$Status),
    Failing = paste0(fmt_num(b$Fail_Percentage, 1), "%"),
    stringsAsFactors = FALSE
  )
}

#' Print cIPMA Results
#'
#' Prints the combined table: importance, performance, the CE-FDH effect size and
#' its permutation p-value, the necessity decision at \code{pi = 1}, the required
#' level at the target level and the share of cases below the required level at
#' the chosen tolerance.
#'
#' @param x A cipma object.
#' @param ... Ignored.
#'
#' @return The object \code{x}, invisibly.
#'
#' @export
print.cipma <- function(x, ...) {
  pi_txt <- pi_label(x$Tolerance)
  cat(sprintf("cIPMA: target %s at level %s (n = %d, tolerance pi = %s)\n\n",
              x$Dependent, format(x$Target_Level), x$Settings$n, pi_txt))
  print(cipma_main_table(x), row.names = FALSE, right = TRUE)
  cat(sprintf(paste0("\nNecessary: CE-FDH d >= 0.10 and p < 0.05 at pi = 1.00. Required and Failing:\n",
                     "bottleneck at pi = %s (NN: no requirement; NA: level not attained).\n"),
              pi_txt))
  cat("summary(x) adds the readings; summary(x, detail = \"full\") adds the tolerance path and\n",
      "the bottleneck table; cipma_diagnostics(x) runs the QNCA diagnostics.\n", sep = "")
  invisible(x)
}

utils::globalVariables(c("Importance", "Performance", "PlotSize", "ColorGroup", "Label",
                         "Fail_Percentage", "Is_Necessary", "Construct", "pi", "d",
                         "Required_Level"))
