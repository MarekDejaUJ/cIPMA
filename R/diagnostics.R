#' QNCA Diagnostics of a cIPMA
#'
#' Runs the diagnostics of Quantile Necessary Condition Analysis on every
#' predictor of a cIPMA, with the rescaled scores and the data-range scope of
#' each pair. The reference band compares the observed effect size with a
#' Gaussian-copula benchmark that keeps the marginals and the rank dependence of
#' the pair and imposes no attainability floor. The deletion screen reports the
#' change of the effect size when one case, or \code{k} cases jointly, leave the
#' sample. The resolution table gives, target by target, the conditioning-set
#' size, the selected rank, the inherited flag and the resistance of the frontier
#' at the chosen tolerance. The consistency probe follows the effect size across
#' subsample sizes. The addition screen, run when \code{points} is given, adds
#' cases at declared positions.
#'
#' @param x A cipma object.
#' @param constructs Predictors to diagnose; default all.
#' @param pi_grid Tolerances of the band, the deletion screen and the addition
#'   screen. Default 1, 0.95, 0.90.
#' @param k Number of cases deleted jointly in the deletion screen. Default 1.
#' @param M Number of benchmark datasets of the reference band. Default 199.
#' @param probe Logical. Run the consistency probe. Default \code{TRUE}.
#' @param reps Subsamples per size of the consistency probe. Default 50.
#' @param points Optional \code{c(x, y)} position, or two-column matrix of
#'   positions, of cases added by the addition screen.
#' @param counts Numbers of added cases of the addition screen. Default 0 to 5.
#' @param n_grid Number of outcome levels. Default 50.
#' @param seed Optional seed; the random number state of the session is
#'   restored afterwards.
#'
#' @return An object of class \code{cipma_diagnostics}: per construct, the
#'   reference band, the deletion screen, the resolution table, the consistency
#'   probe and the addition screen, with the settings.
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
#' diag <- cipma_diagnostics(res, M = 49, reps = 10, seed = 1)
#' diag
#' }
#' @export
cipma_diagnostics <- function(x, constructs = NULL, pi_grid = c(1, 0.95, 0.90), k = 1L,
                              M = 199L, probe = TRUE, reps = 50L, points = NULL,
                              counts = 0:5, n_grid = 50L, seed = NULL) {
  if (!inherits(x, "cipma")) stop("'x' must be a cipma object.")
  if (is.null(constructs)) constructs <- x$Predictors
  unknown <- setdiff(constructs, x$Predictors)
  if (length(unknown)) stop("Unknown constructs: ", paste(unknown, collapse = ", "))
  Y <- x$Data_Rescaled[[x$Dependent]]
  res <- list()
  with_seed(seed, {
    for (i in seq_along(constructs)) {
      cn <- constructs[i]
      X <- x$Data_Rescaled[[cn]]
      ok <- is.finite(X) & is.finite(Y)
      xx <- X[ok]
      yy <- Y[ok]
      s <- if (is.null(seed)) NULL else seed + i
      one <- list(
        band = qnca::spuriousness_band(xx, yy, pi_grid = pi_grid, n_grid = n_grid, M = M, seed = s)$band,
        outliers = qnca::qnca_outliers(xx, yy, pi_grid = pi_grid, k = k, n_grid = n_grid),
        resolution = qnca::qnca_resolution(xx, yy, pi = x$Tolerance, n_grid = n_grid)
      )
      if (isTRUE(probe)) {
        one$probe <- qnca::consistency_probe(xx, yy, pi_pair = c(1, min(pi_grid)), reps = reps,
                                             n_grid = n_grid, seed = s)
      }
      if (!is.null(points)) {
        one$sensitivity <- qnca::qnca_sensitivity(xx, yy, points = points, counts = counts,
                                                  pi_grid = pi_grid, n_grid = n_grid)
      }
      res[[cn]] <- one
    }
  })
  out <- list(results = res, dependent = x$Dependent, target = x$Target_Level,
              tolerance = x$Tolerance,
              settings = list(pi_grid = pi_grid, k = k, M = M, probe = probe, reps = reps,
                              points = points, counts = counts, n_grid = n_grid, seed = seed))
  class(out) <- "cipma_diagnostics"
  out
}

# Tables of the diagnostics, one row per construct and tolerance.
diagnostics_tables <- function(x) {
  band <- do.call(rbind, lapply(names(x$results), function(cn) {
    b <- x$results[[cn]]$band
    data.frame(Construct = cn, pi = pi_label(b$pi), d_obs = fmt_num(b$d_obs), upper = fmt_num(b$upper),
               excess = fmt_num(b$excess), p_null = fmt_num(b$p_null),
               Above_band = fmt_yes(b$excess > 0), stringsAsFactors = FALSE)
  }))
  screen <- do.call(rbind, lapply(names(x$results), function(cn) {
    o <- x$results[[cn]]$outliers$table
    rel <- grep("^dif_rel_", names(o), value = TRUE)
    top <- o[1, , drop = FALSE]
    data.frame(Construct = cn, Cases = top$cases,
               Relative_change = paste(vapply(rel, function(r) fmt_num(top[[r]]), ""), collapse = " / "),
               Flagged_rows = sum(o$flagged), stringsAsFactors = FALSE)
  }))
  resolution <- do.call(rbind, lapply(names(x$results), function(cn) {
    r <- x$results[[cn]]$resolution
    upto <- r[r$y <= x$target, , drop = FALSE]
    data.frame(Construct = cn, Targets = nrow(upto),
               Inherited = sum(upto$inherited),
               Smallest_resistance = if (nrow(upto)) min(upto$resistance) else NA_integer_,
               Smallest_attaining_set = if (nrow(upto)) min(upto$k) else NA_integer_,
               stringsAsFactors = FALSE)
  }))
  probe <- NULL
  if (isTRUE(x$settings$probe)) {
    probe <- do.call(rbind, lapply(names(x$results), function(cn) {
      p <- x$results[[cn]]$probe
      data.frame(Construct = cn, Divergence = fmt_num(p$divergence),
                 Drift = paste(fmt_num(p$beta_drift, 5), collapse = " / "), stringsAsFactors = FALSE)
    }))
  }
  list(band = band, screen = screen, resolution = resolution, probe = probe)
}

#' @export
print.cipma_diagnostics <- function(x, ...) {
  t <- diagnostics_tables(x)
  s <- x$settings
  cat(sprintf("QNCA diagnostics of the cIPMA predictors of %s (tolerances %s)\n\n",
              x$dependent, paste(pi_label(s$pi_grid), collapse = ", ")))
  cat(sprintf("Reference band (%d benchmark datasets): observed d, upper band limit, excess\n", s$M))
  print(t$band, row.names = FALSE, right = TRUE)
  cat(sprintf("\nDeletion screen (k = %d): most influential deletion and its relative change of d per tolerance\n", s$k))
  print(t$screen, row.names = FALSE, right = TRUE)
  cat(sprintf("\nResolution at pi = %s for outcome targets up to %s\n", pi_label(x$tolerance), format(x$target)))
  print(t$resolution, row.names = FALSE, right = TRUE)
  if (!is.null(t$probe)) {
    cat("\nConsistency probe: divergence of the drift slopes and the slopes at the two tolerances\n")
    print(t$probe, row.names = FALSE, right = TRUE)
  }
  cat("\nAn empty zone above the band exceeds the dependence-matched benchmark; the band is a\n",
      "proposed diagnostic with operating characteristics still to be established.\n", sep = "")
  invisible(x)
}

diagnostics_markdown <- function(x, md_table) {
  t <- diagnostics_tables(x)
  s <- x$settings
  out <- c("## QNCA diagnostics", "",
           sprintf("Reference band with %d benchmark datasets:", s$M), "", md_table(t$band), "",
           sprintf("Deletion screen with k = %d:", s$k), "", md_table(t$screen), "",
           sprintf("Resolution at pi = %s for outcome targets up to %s:", pi_label(x$tolerance), format(x$target)),
           "", md_table(t$resolution))
  if (!is.null(t$probe)) out <- c(out, "", "Consistency probe:", "", md_table(t$probe))
  out
}
