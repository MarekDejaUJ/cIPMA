#' Summarise cIPMA Results
#'
#' Collects the results of a cIPMA at three levels of detail. \code{"brief"}
#' gives the combined table. \code{"standard"} adds the PLS-SEM p-values, the
#' effect size and p-value at the chosen tolerance and one reading per
#' construct. \code{"full"} adds the tolerance path of the effect size and of the
#' required level, the bottleneck table across outcome levels and the settings.
#'
#' @param object A cipma object.
#' @param detail One of \code{"standard"}, \code{"brief"} or \code{"full"}.
#' @param digits Decimals of effect sizes and importance.
#' @param ... Ignored.
#'
#' @return An object of class \code{summary.cipma}, printed by its print method.
#'
#' @examples
#' \donttest{
#' library(seminr)
#' mm <- constructs(
#'   composite("Image", multi_items("IMAG", 1:5)),
#'   composite("Value", multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3))
#' )
#' sm <- relationships(paths(from = c("Image", "Value"), to = "Satisfaction"))
#' pls <- estimate_pls(mobi, mm, sm)
#' boot <- bootstrap_model(pls, nboot = 50, cores = 1, seed = 1)
#' scales <- setNames(rep(list(c(1, 10)), ncol(mobi)), names(mobi))
#' res <- cipma(boot, "Satisfaction", mobi, scales, target_level = 80,
#'              pi = 0.95, B = 99, seed = 1)
#' summary(res, detail = "full")
#' }
#' @export
summary.cipma <- function(object, detail = c("standard", "brief", "full"), digits = 3, ...) {
  detail <- match.arg(detail)
  x <- object
  out <- list(detail = detail, header = list(
    dependent = x$Dependent, target = x$Target_Level, pi = x$Tolerance,
    n = x$Settings$n, B = x$Settings$B, n_grid = x$Settings$n_grid,
    n_grid_test = x$Settings$n_grid_test, seed = x$Settings$seed,
    qnca_version = x$Settings$qnca_version))
  out$main <- cipma_main_table(x, digits)
  r <- x$cIPMA_results
  b <- x$Bottleneck_summary[match(r$Construct, x$Bottleneck_summary$Construct), ]
  if (detail %in% c("standard", "full")) {
    out$detail_table <- data.frame(
      Construct = r$Construct,
      PLS_p = fmt_p(r$PLS_p_value),
      d_pi = fmt_num(r$d_pi, digits),
      p_pi = fmt_p(r$p_pi),
      Failing_cases = format(b$Fail_Count),
      stringsAsFactors = FALSE)
    names(out$detail_table)[3:4] <- paste0(c("d_", "p_"), pi_label(x$Tolerance))
    out$readings <- cipma_readings(x)
  }
  if (detail == "full") {
    path <- x$Tolerance_Path
    out$path_d <- path_wide(path, "d", function(v) fmt_num(v, digits))
    path$Required_fmt <- fmt_required(path$Required_Level, path$Status)
    out$path_required <- path_wide(path, "Required_fmt", identity)
    full <- x$Bottleneck_full
    full$Required_fmt <- fmt_required(full$Required_Level, full$Status)
    out$bottleneck <- level_wide(full)
  }
  class(out) <- "summary.cipma"
  out
}

#' @export
print.summary.cipma <- function(x, ...) {
  h <- x$header
  cat(sprintf("Combined importance-performance map analysis\n"))
  cat(sprintf("Target construct: %s, target level: %s, cases: %d\n", h$dependent, format(h$target), h$n))
  cat(sprintf("Tolerance of the bottleneck: pi = %s; permutations: %s; outcome levels: %d (test: %d); qnca %s\n\n",
              pi_label(h$pi), format(h$B), h$n_grid, h$n_grid_test, h$qnca_version))
  cat("Combined map\n")
  print(x$main, row.names = FALSE, right = TRUE)
  if (!is.null(x$detail_table)) {
    cat("\nPLS-SEM significance, effect at the chosen tolerance, failing cases\n")
    print(x$detail_table, row.names = FALSE, right = TRUE)
    cat("\nReadings\n")
    cat(paste0("- ", x$readings), sep = "\n")
    if (h$pi < 1) {
      cat(sprintf(paste0("\nEffect sizes at pi = %s are relative to the tolerance and rise as pi falls;\n",
                         "the NCA landmark of 0.10 applies to the CE-FDH effect size at pi = 1.00.\n"),
                  pi_label(h$pi)))
    }
  }
  if (!is.null(x$path_d)) {
    cat("\nTolerance path: effect size d\n")
    print(x$path_d, row.names = FALSE, right = TRUE)
    cat(sprintf("\nTolerance path: required level for %s >= %s\n", h$dependent, format(h$target)))
    print(x$path_required, row.names = FALSE, right = TRUE)
    cat(sprintf("\nBottleneck table at pi = %s (NN: no requirement; NA: level not attained)\n", pi_label(h$pi)))
    print(x$bottleneck, row.names = FALSE, right = TRUE)
  }
  invisible(x)
}

# One sentence per construct in the voice of the guide.
cipma_readings <- function(x) {
  r <- x$cIPMA_results
  b <- x$Bottleneck_summary[match(r$Construct, x$Bottleneck_summary$Construct), ]
  vapply(seq_len(nrow(r)), function(i) {
    nec <- if (is.na(r$Is_Necessary[i])) "carries no necessity decision (no permutation test)"
      else if (r$Is_Necessary[i]) "is necessary at the strict limit"
      else "is not necessary at the strict limit"
    req <- switch(b$Status[i],
      "required" = sprintf("the required level of %s for %s >= %s at pi = %s is %s, with %s%% of the cases below the required level",
                           r$Construct[i], x$Dependent, format(x$Target_Level), pi_label(x$Tolerance),
                           fmt_num(b$Required_Level[i], 2), fmt_num(b$Fail_Percentage[i], 1)),
      "not necessary" = sprintf("the target level of %s sets no requirement on %s at pi = %s",
                                x$Dependent, r$Construct[i], pi_label(x$Tolerance)),
      "not attained" = sprintf("no case attains the target level of %s", x$Dependent))
    sprintf("%s %s (d = %s, p = %s); %s.", r$Construct[i], nec, fmt_num(r$d[i], 3), fmt_p(r$p_value[i]), req)
  }, character(1))
}

path_wide <- function(path, column, fmt) {
  pis <- sort(unique(path$pi), decreasing = TRUE)
  out <- data.frame(pi = pi_label(pis), stringsAsFactors = FALSE)
  for (cn in unique(path$Construct)) {
    s <- path[path$Construct == cn, ]
    out[[cn]] <- fmt(s[[column]][match(pis, s$pi)])
  }
  out
}

level_wide <- function(full) {
  lv <- sort(unique(full$Outcome_Level))
  out <- data.frame(Outcome_Level = format(lv), stringsAsFactors = FALSE)
  for (cn in unique(full$Construct)) {
    s <- full[full$Construct == cn, ]
    out[[cn]] <- s$Required_fmt[match(lv, s$Outcome_Level)]
  }
  out
}

#' Write a Full cIPMA Report
#'
#' Writes the full summary of a cIPMA, and of its diagnostics when given, as a
#' Markdown report.
#'
#' @param x A cipma object.
#' @param file Optional path of the Markdown file; the report is returned only
#'   when \code{file} is \code{NULL}.
#' @param diagnostics Optional result of \code{\link{cipma_diagnostics}}.
#'
#' @return The lines of the report, invisibly.
#'
#' @examples
#' \donttest{
#' library(seminr)
#' mm <- constructs(
#'   composite("Image", multi_items("IMAG", 1:5)),
#'   composite("Value", multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3))
#' )
#' sm <- relationships(paths(from = c("Image", "Value"), to = "Satisfaction"))
#' pls <- estimate_pls(mobi, mm, sm)
#' boot <- bootstrap_model(pls, nboot = 50, cores = 1, seed = 1)
#' scales <- setNames(rep(list(c(1, 10)), ncol(mobi)), names(mobi))
#' res <- cipma(boot, "Satisfaction", mobi, scales, target_level = 80, B = 99, seed = 1)
#' report <- cipma_report(res, file = file.path(tempdir(), "cipma-report.md"))
#' }
#' @export
cipma_report <- function(x, file = NULL, diagnostics = NULL) {
  s <- summary(x, detail = "full")
  h <- s$header
  md_table <- function(df) {
    c(paste0("| ", paste(names(df), collapse = " | "), " |"),
      paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
      apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  }
  lines <- c(
    sprintf("# cIPMA report: %s at level %s", h$dependent, format(h$target)), "",
    sprintf("Cases: %d. Tolerance of the bottleneck: pi = %s. Permutations: %s. Outcome levels: %d (test: %d). qnca %s.",
            h$n, pi_label(h$pi), format(h$B), h$n_grid, h$n_grid_test, h$qnca_version), "",
    "## Combined map", "", md_table(s$main), "",
    "Necessary: CE-FDH d >= 0.10 and p < 0.05 at pi = 1.00. Required and Failing: bottleneck at the chosen tolerance.", "",
    "## PLS-SEM significance and effect at the chosen tolerance", "", md_table(s$detail_table), "",
    "## Readings", "", paste0("- ", s$readings), "",
    "## Tolerance path: effect size d", "", md_table(s$path_d), "",
    sprintf("## Tolerance path: required level for %s >= %s", h$dependent, format(h$target)), "",
    md_table(s$path_required), "",
    sprintf("## Bottleneck table at pi = %s", pi_label(h$pi)), "", md_table(s$bottleneck), "",
    "NN: no requirement. NA: level not attained.")
  if (!is.null(diagnostics)) lines <- c(lines, "", diagnostics_markdown(diagnostics, md_table))
  if (!is.null(file)) writeLines(lines, file, useBytes = TRUE)
  invisible(lines)
}
