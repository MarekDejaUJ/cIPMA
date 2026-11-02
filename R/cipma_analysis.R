#' Run Combined Importance-Performance Map Analysis (cIPMA)
#'
#' Combines the importance-performance map analysis of a PLS-SEM model with
#' Quantile Necessary Condition Analysis (QNCA) of each predictor of the target
#' construct (Hauff et al. 2024). Importance is the total effect on the target
#' construct, and performance is the mean of the latent variable scores rescaled
#' to 0-100. Necessity is decided at the strict limit \code{pi = 1}, where the
#' QNCA frontier equals the CE-FDH ceiling of Necessary Condition Analysis: a
#' predictor is necessary when its CE-FDH effect size reaches 0.10 and its
#' permutation p-value stays below 0.05. The tolerance \code{pi} sets the share
#' \code{1 - pi} of cases that attain the target level and may sit below the
#' reported requirement; the bottleneck, the share of cases failing it and the
#' effect size at \code{pi} follow the chosen tolerance. The tolerance path
#' repeats the effect size and the bottleneck for every tolerance in
#' \code{pi_path}.
#'
#' @param model A bootstrapped seminr model object.
#' @param target_construct String. The name of the target construct.
#' @param data The raw data frame used for the model.
#' @param scales A named list of theoretical minimum and maximum per indicator.
#' @param target_level Numeric. The desired level of the target construct on the
#'   0-100 scale. Default 85.
#' @param nca_rep Deprecated. Number of permutations; use \code{B}.
#' @param pls_model Optional. The base (non-bootstrapped) seminr model for the
#'   assumption checks.
#' @param pi Tolerance in (0, 1] for the bottleneck and the effect size at the
#'   chosen tolerance. Default 1, which reproduces the CE-FDH results of NCA.
#' @param pi_path Tolerances of the tolerance path. Default 1.00 to 0.90 in steps
#'   of 0.01.
#' @param B Number of permutations of the significance test. Default 1999; 0
#'   skips the test.
#' @param n_grid Number of outcome levels of the QNCA effect sizes. Default 400.
#' @param n_grid_test Number of outcome levels inside the permutation test.
#'   Default 100.
#' @param seed Optional seed of the permutation test; the random number state of
#'   the session is restored afterwards.
#'
#' @return An object of class \code{cipma}: a list with the combined results
#'   (\code{cIPMA_results}), the bottleneck at the target level
#'   (\code{Bottleneck_summary}), the bottleneck table across outcome levels
#'   (\code{Bottleneck_full}), the tolerance path (\code{Tolerance_Path}), the
#'   settings, the rescaled scores and the models.
#'
#' @references Hauff, S., Richter, N. F., Sarstedt, M., & Ringle, C. M. (2024).
#'   Importance and performance in PLS-SEM and NCA: Introducing the combined
#'   importance-performance map analysis (cIPMA). Journal of Retailing and
#'   Consumer Services, 78, 103723. \doi{10.1016/j.jretconser.2024.103723}
#'
#' @examples
#' \donttest{
#' library(seminr)
#' mm <- constructs(
#'   composite("Image", multi_items("IMAG", 1:5)),
#'   composite("Expectation", multi_items("CUEX", 1:3)),
#'   composite("Value", multi_items("PERV", 1:2)),
#'   composite("Satisfaction", multi_items("CUSA", 1:3))
#' )
#' sm <- relationships(
#'   paths(from = c("Image", "Expectation", "Value"), to = "Satisfaction")
#' )
#' pls <- estimate_pls(mobi, mm, sm)
#' boot <- bootstrap_model(pls, nboot = 50, cores = 1, seed = 1)
#' scales <- setNames(rep(list(c(1, 10)), ncol(mobi)), names(mobi))
#' res <- cipma(boot, "Satisfaction", mobi, scales, target_level = 80,
#'              pi = 0.95, B = 99, seed = 1, pls_model = pls)
#' res
#' summary(res)
#' }
#' @importFrom stats pnorm sd
#' @export
cipma <- function(model, target_construct, data, scales, target_level = 85,
                  nca_rep = NULL, pls_model = NULL, pi = 1,
                  pi_path = seq(1, 0.90, by = -0.01), B = 1999L, n_grid = 400L,
                  n_grid_test = 100L, seed = NULL) {

  if (is.null(model$boot_paths)) {
    stop("The seminr model must be bootstrapped to obtain CIs and p-values.")
  }
  if (!is.null(nca_rep)) {
    warning("'nca_rep' is deprecated; the number of permutations is now set by 'B'.")
    B <- nca_rep
  }
  if (length(pi) != 1L || !is.finite(pi) || pi <= 0 || pi > 1) {
    stop("'pi' must be one number in (0, 1].")
  }
  if (length(pi_path) < 1L || any(!is.finite(pi_path)) || any(pi_path <= 0 | pi_path > 1)) {
    stop("'pi_path' must hold numbers in (0, 1].")
  }
  if (length(target_level) != 1L || !is.finite(target_level) ||
      target_level < 0 || target_level > 100) {
    stop("'target_level' must be one number in [0, 100].")
  }
  pi_path <- sort(unique(round(c(pi_path, 1, pi), 9)), decreasing = TRUE)

  smry <- summary(model)
  boot_total <- smry$bootstrapped_total_paths

  pattern <- paste0("\\s*->\\s*", target_construct, "$")
  target_rows_idx <- grep(pattern, rownames(boot_total))

  if (length(target_rows_idx) == 0) {
    stop(paste0("Target construct '", target_construct, "' not found as an outcome."))
  }

  target_matrix <- boot_total[target_rows_idx, , drop = FALSE]
  predictors <- sub(pattern, "", rownames(target_matrix))

  pls_res <- data.frame(
    Construct = predictors,
    Importance = NA,
    Performance = NA,
    PLS_p_value = NA,
    stringsAsFactors = FALSE
  )

  original_indicators <- as.data.frame(data)
  std_weights_matrix <- model$outer_weights
  std_weights_df <- as.data.frame(as.table(std_weights_matrix))
  colnames(std_weights_df) <- c("Indicator", "Latent_variable", "Standardized_weight")
  std_weights_df <- std_weights_df[std_weights_df$Standardized_weight != 0, ]

  negative_weights <- std_weights_df$Standardized_weight < 0
  if (any(negative_weights)) {
    neg_indicators <- std_weights_df$Indicator[negative_weights]
    warning(paste0(
      "Negative outer weights detected for: ",
      paste(neg_indicators, collapse = ", "),
      ". This may invert performance scores. Consider checking for reverse-coded items."
    ))
  }

  if (!is.null(model$sdData)) {
    std_weights_df$SD <- model$sdData[as.character(std_weights_df$Indicator)]
  } else {
    std_weights_df$SD <- apply(data[, as.character(std_weights_df$Indicator)], 2, stats::sd, na.rm = TRUE)
  }

  std_weights_df$Unstandardized_weight <- std_weights_df$Standardized_weight / std_weights_df$SD

  sum_weights <- tapply(std_weights_df$Unstandardized_weight, std_weights_df$Latent_variable, sum)
  std_weights_df$Sum_weight <- sum_weights[as.character(std_weights_df$Latent_variable)]
  std_weights_df$Norm_Unstd_weight <- std_weights_df$Unstandardized_weight / std_weights_df$Sum_weight

  lv_unscaled <- data.frame(row.names = 1:nrow(data))
  unique_lvs <- unique(std_weights_df$Latent_variable)

  for (lv in unique_lvs) {
    lv_dat <- std_weights_df[std_weights_df$Latent_variable == lv, ]
    inds <- as.character(lv_dat$Indicator)
    wts <- lv_dat$Norm_Unstd_weight
    ind_scores <- as.matrix(original_indicators[, inds, drop = FALSE])
    lv_unscaled[[as.character(lv)]] <- as.vector(ind_scores %*% wts)
  }

  lv_df <- lv_unscaled
  for (lv in names(lv_unscaled)) {
    related_ind <- as.character(std_weights_df$Indicator[std_weights_df$Latent_variable == lv][1])
    if (related_ind %in% names(scales)) {
      min_val <- scales[[related_ind]][1]
      max_val <- scales[[related_ind]][2]
      norm_val <- ((lv_unscaled[[lv]] - min_val) / (max_val - min_val)) * 100
      norm_val <- pmax(0, pmin(100, norm_val))
      lv_df[[lv]] <- norm_val
    }
  }

  for (i in 1:nrow(pls_res)) {
    pred_name <- pls_res$Construct[i]
    row_key <- rownames(target_matrix)[grep(paste0("^", pred_name, "\\s*->"), rownames(target_matrix))][1]

    pls_res$Importance[i] <- target_matrix[row_key, "Original Est."]

    cols <- colnames(target_matrix)
    if ("p" %in% cols) {
      pls_res$PLS_p_value[i] <- target_matrix[row_key, "p"]
    } else if ("T Stat." %in% cols) {
      pls_res$PLS_p_value[i] <- 2 * (1 - stats::pnorm(abs(target_matrix[row_key, "T Stat."])))
    }

    if (pred_name %in% colnames(lv_df)) {
      pls_res$Performance[i] <- mean(lv_df[[pred_name]], na.rm = TRUE)
    }
  }

  Y <- lv_df[[target_construct]]
  levels_full <- sort(unique(c(seq(0, 100, by = 10), target_level)))

  nec <- data.frame(Construct = predictors, d = NA_real_, p_value = NA_real_,
                    Is_Necessary = NA, d_pi = NA_real_, p_pi = NA_real_,
                    stringsAsFactors = FALSE)
  bn <- data.frame(Construct = predictors, Required_Level = NA_real_, Status = NA_character_,
                   Fail_Count = NA_real_, Fail_Percentage = NA_real_, stringsAsFactors = FALSE)
  path_rows <- list()
  full_rows <- list()

  with_seed(seed, {
    for (i in seq_along(predictors)) {
      pred <- predictors[i]
      X <- lv_df[[pred]]
      ok <- is.finite(X) & is.finite(Y)
      x <- X[ok]
      y <- Y[ok]

      nec$d[i] <- qnca::nca_ce_fdh_d(x, y)
      nec$p_value[i] <- qnca_pvalue(x, y, 1, n_grid_test, B)
      if (!is.na(nec$p_value[i])) {
        nec$Is_Necessary[i] <- nec$d[i] >= 0.10 && nec$p_value[i] < 0.05
      }
      if (pi < 1) {
        nec$d_pi[i] <- qnca_effect(x, y, pi, n_grid)
        nec$p_pi[i] <- qnca_pvalue(x, y, pi, n_grid_test, B)
      } else {
        nec$d_pi[i] <- nec$d[i]
        nec$p_pi[i] <- nec$p_value[i]
      }

      req <- qnca_requirement(x, y, pi, target_level)
      fail <- qnca_failing(x, req$required)
      bn$Required_Level[i] <- req$required
      bn$Status[i] <- req$status
      bn$Fail_Count[i] <- fail[["count"]]
      bn$Fail_Percentage[i] <- fail[["percent"]]

      for (p in pi_path) {
        rq <- qnca_requirement(x, y, p, target_level)
        fl <- qnca_failing(x, rq$required)
        path_rows[[length(path_rows) + 1L]] <- data.frame(
          Construct = pred, pi = p,
          d = if (p == 1) nec$d[i] else qnca_effect(x, y, p, n_grid),
          Required_Level = rq$required, Status = rq$status,
          Fail_Percentage = fl[["percent"]], stringsAsFactors = FALSE)
      }

      for (lvl in levels_full) {
        rq <- qnca_requirement(x, y, pi, lvl)
        full_rows[[length(full_rows) + 1L]] <- data.frame(
          Construct = pred, Outcome_Level = lvl, Required_Level = rq$required,
          Status = rq$status, stringsAsFactors = FALSE)
      }
    }
  })

  results <- merge(pls_res, nec, by = "Construct", sort = FALSE)
  results <- merge(results, bn[, c("Construct", "Required_Level", "Fail_Percentage")],
                   by = "Construct", sort = FALSE)
  results <- results[match(predictors, results$Construct), ]
  rownames(results) <- NULL

  output <- list(
    cIPMA_results = results,
    Bottleneck_summary = bn,
    Bottleneck_full = do.call(rbind, full_rows),
    Tolerance_Path = do.call(rbind, path_rows),
    Target_Level = target_level,
    Tolerance = pi,
    Settings = list(pi = pi, pi_path = pi_path, B = B, n_grid = n_grid,
                    n_grid_test = n_grid_test, seed = seed,
                    n = sum(is.finite(Y)),
                    qnca_version = as.character(utils::packageVersion("qnca"))),
    Predictors = predictors,
    Dependent = target_construct,
    Data_Rescaled = lv_df,
    Boot_Model = model,
    PLS_Model = pls_model
  )

  class(output) <- "cipma"
  return(output)
}

#' Check cIPMA Assumptions and Model Diagnostics
#'
#' @param x A cipma object or a seminr model object.
#' @param vif_threshold Numeric. VIF threshold for multicollinearity. Default is 5.
#' @param htmt_threshold Numeric. HTMT threshold for discriminant validity. Default is 0.90.
#'
#' @return A list containing assumption checks and diagnostics.
#'
#' @export
check_assumptions <- function(x, vif_threshold = 5, htmt_threshold = 0.90) {

  if (inherits(x, "cipma")) {
    if (!is.null(x$PLS_Model)) {
      model <- x$PLS_Model
    } else {
      model <- x$Boot_Model
    }
  } else if (inherits(x, "boot_seminr_model") || inherits(x, "seminr_model")) {
    model <- x
  } else {
    stop("Input must be a cipma object or a seminr model.")
  }

  model_smry <- summary(model)

  outer_weights <- model$outer_weights
  weights_df <- as.data.frame(as.table(outer_weights))
  colnames(weights_df) <- c("Indicator", "Construct", "Weight")
  weights_df <- weights_df[weights_df$Weight != 0, ]

  negative_weights <- weights_df[weights_df$Weight < 0, ]
  weights_check <- list(
    passed = nrow(negative_weights) == 0,
    negative_indicators = if (nrow(negative_weights) > 0) negative_weights else NULL
  )

  reliability <- tryCatch({
    model_smry$reliability
  }, error = function(e) NULL)

  if (!is.null(reliability)) {
    alpha_check <- reliability[, "alpha"] >= 0.70
    rhoC_check <- reliability[, "rhoC"] >= 0.70
    ave_check <- reliability[, "AVE"] >= 0.50

    problematic <- character(0)
    for (i in seq_len(nrow(reliability))) {
      if (is.na(alpha_check[i]) || is.na(rhoC_check[i]) || is.na(ave_check[i])) next
      if (!alpha_check[i] || !rhoC_check[i] || !ave_check[i]) {
        problematic <- c(problematic, rownames(reliability)[i])
      }
    }

    reliability_check <- list(
      table = reliability,
      alpha_passed = all(alpha_check, na.rm = TRUE),
      rhoC_passed = all(rhoC_check, na.rm = TRUE),
      ave_passed = all(ave_check, na.rm = TRUE),
      problematic_constructs = problematic
    )
  } else {
    reliability_check <- list(
      table = NULL,
      note = "Reliability not available."
    )
  }

  htmt <- tryCatch({
    model_smry$validity$htmt
  }, error = function(e) NULL)

  if (!is.null(htmt)) {
    htmt_numeric <- htmt
    htmt_numeric[htmt_numeric == "."] <- NA
    htmt_numeric <- apply(htmt_numeric, 2, as.numeric)
    rownames(htmt_numeric) <- rownames(htmt)

    htmt_values <- htmt_numeric[lower.tri(htmt_numeric)]
    htmt_exceeded <- any(htmt_values > htmt_threshold, na.rm = TRUE)

    problematic_pairs <- character(0)
    if (htmt_exceeded && nrow(htmt_numeric) > 1) {
      for (i in 2:nrow(htmt_numeric)) {
        for (j in 1:(i-1)) {
          val <- htmt_numeric[i, j]
          if (!is.na(val) && val > htmt_threshold) {
            problematic_pairs <- c(problematic_pairs,
                                   paste(rownames(htmt_numeric)[i], "<->", colnames(htmt_numeric)[j]))
          }
        }
      }
    }

    htmt_check <- list(
      table = htmt,
      passed = !htmt_exceeded,
      threshold = htmt_threshold,
      problematic_pairs = if (length(problematic_pairs) > 0) problematic_pairs else NULL
    )
  } else {
    htmt_check <- list(
      table = NULL,
      passed = NA,
      threshold = htmt_threshold,
      note = "HTMT not available."
    )
  }

  fl_criteria <- tryCatch({
    model_smry$validity$fl_criteria
  }, error = function(e) NULL)

  vif_antecedents <- tryCatch({
    model_smry$vif_antecedents
  }, error = function(e) NULL)

  if (!is.null(vif_antecedents)) {
    all_vifs <- unlist(vif_antecedents)
    vif_exceeded <- any(all_vifs > vif_threshold, na.rm = TRUE)
    problematic_vif <- names(all_vifs)[!is.na(all_vifs) & all_vifs > vif_threshold]

    vif_check <- list(
      table = vif_antecedents,
      passed = !vif_exceeded,
      threshold = vif_threshold,
      problematic = if (length(problematic_vif) > 0) problematic_vif else NULL
    )
  } else {
    vif_check <- list(
      table = NULL,
      passed = NA,
      threshold = vif_threshold,
      note = "VIF not available."
    )
  }

  sample_size <- nrow(model$data)
  n_constructs <- ncol(model$construct_scores)
  min_recommended <- max(50, n_constructs * 10)

  sample_check <- list(
    n = sample_size,
    min_recommended = min_recommended,
    passed = sample_size >= min_recommended
  )

  output <- list(
    positive_weights = weights_check,
    reliability = reliability_check,
    discriminant_validity_htmt = htmt_check,
    fornell_larcker = fl_criteria,
    multicollinearity_vif = vif_check,
    sample_size = sample_check
  )

  class(output) <- "cipma_assumptions"
  return(output)
}



#' Print Assumption Check Results
#'
#' @param x A cipma_assumptions object.
#' @param ... Ignored.
#'
#' @export
print.cipma_assumptions <- function(x, ...) {

  cat("=== cIPMA Assumption Checks ===\n\n")

  cat("1. POSITIVE OUTER WEIGHTS\n")
  if (x$positive_weights$passed) {
    cat("   [PASS] All outer weights are positive.\n")
  } else {
    cat("   [FAIL] Negative weights detected:\n")
    print(x$positive_weights$negative_indicators)
  }

  cat("\n2. RELIABILITY (alpha >= 0.70, rhoC >= 0.70, AVE >= 0.50)\n")
  if (!is.null(x$reliability$table)) {
    print(round(x$reliability$table, 3))
    if (length(x$reliability$problematic_constructs) > 0) {
      cat("   [WARN] Check constructs:", paste(x$reliability$problematic_constructs, collapse = ", "), "\n")
    } else {
      cat("   [PASS] All reliability thresholds met.\n")
    }
  } else {
    cat("   [INFO]", x$reliability$note, "\n")
  }

  cat("\n3. DISCRIMINANT VALIDITY (HTMT <", x$discriminant_validity_htmt$threshold, ")\n")
  if (!is.null(x$discriminant_validity_htmt$table)) {
    print(x$discriminant_validity_htmt$table, quote = FALSE)
    if (x$discriminant_validity_htmt$passed) {
      cat("   [PASS] All HTMT values below threshold.\n")
    } else {
      cat("   [FAIL] HTMT exceeded for:\n")
      for (pair in x$discriminant_validity_htmt$problematic_pairs) {
        cat("         ", pair, "\n")
      }
    }
  } else {
    cat("   [INFO]", x$discriminant_validity_htmt$note, "\n")
  }

  cat("\n4. MULTICOLLINEARITY (VIF <", x$multicollinearity_vif$threshold, ")\n")
  if (!is.null(x$multicollinearity_vif$table)) {
    for (dv in names(x$multicollinearity_vif$table)) {
      cat("  ", dv, ":\n")
      vifs <- x$multicollinearity_vif$table[[dv]]
      cat("   ", paste(names(vifs), "=", round(vifs, 2), collapse = ", "), "\n")
    }
    if (x$multicollinearity_vif$passed) {
      cat("   [PASS] No multicollinearity issues detected.\n")
    } else {
      cat("   [FAIL] High VIF for:", paste(x$multicollinearity_vif$problematic, collapse = ", "), "\n")
    }
  } else {
    cat("   [INFO]", x$multicollinearity_vif$note, "\n")
  }

  cat("\n5. SAMPLE SIZE\n")
  cat("   N =", x$sample_size$n, "(recommended >=", x$sample_size$min_recommended, ")\n")
  if (x$sample_size$passed) {
    cat("   [PASS] Sample size adequate.\n")
  } else {
    cat("   [WARN] Sample size may be insufficient.\n")
  }

  cat("\n")
}
