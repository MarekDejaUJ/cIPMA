# Internal helpers: QNCA quantities of one predictor-outcome pair.

# Evaluate `code` with a seed, then restore the random number state of the session.
with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env, inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = env)
    } else if (exists(".Random.seed", envir = env, inherits = FALSE)) {
      rm(".Random.seed", envir = env)
    }
  })
  set.seed(seed)
  force(code)
}

# Required level of the predictor for `level` of the outcome on the QNCA frontier
# at tolerance `pi`. The frontier is evaluated at every observed outcome value up
# to the level, so the monotone envelope at the level is exact. Status follows the
# bottleneck table of NCA: "not necessary" when the requirement equals the lowest
# observed predictor value, "not attained" when no case reaches the level.
qnca_requirement <- function(X, Y, pi, level) {
  if (!any(Y >= level)) {
    return(list(required = NA_real_, status = "not attained"))
  }
  grid <- sort(unique(c(Y[Y <= level], level)))
  phi <- qnca::qnca_frontier(X, Y, pi = pi, y_grid = grid)
  req <- phi[length(phi)]
  if (!is.finite(req) || req <= min(X) + sqrt(.Machine$double.eps)) {
    return(list(required = NA_real_, status = "not necessary"))
  }
  list(required = req, status = "required")
}

# Share of cases below a required level (percent) and their number.
qnca_failing <- function(X, required) {
  if (is.na(required)) return(c(count = 0, percent = 0))
  count <- sum(X < required - sqrt(.Machine$double.eps), na.rm = TRUE)
  c(count = count, percent = 100 * count / sum(!is.na(X)))
}

# Effect size of the QNCA frontier at tolerance `pi` on the data-range scope.
qnca_effect <- function(X, Y, pi, n_grid) {
  qnca::qnca(X, Y, pi = pi, n_grid = n_grid, B = 0L)$d_pi
}

# Permutation p-value of the QNCA effect size at tolerance `pi`.
qnca_pvalue <- function(X, Y, pi, n_grid, B) {
  if (is.na(B) || B < 1) return(NA_real_)
  qnca::qnca(X, Y, pi = pi, n_grid = n_grid, B = as.integer(B))$p_pi
}

# Tolerances written with two decimals for labels and tables.
pi_label <- function(pi) formatC(pi, format = "f", digits = 2)
