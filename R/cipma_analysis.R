#' Run Combined Importance-Performance Map Analysis (cIPMA)
#'
#' @param model A bootstrapped seminr model object.
#' @param target_construct String. The name of the dependent variable.
#' @param data The raw data frame used for the model.
#' @param scales A named list of theoretical min/max for indicators.
#' @param target_level Numeric. The desired level of the target construct (0-100 scale). Default is 85.
#' @param nca_rep Integer. Number of permutations for NCA. Default 10000.
#' @param pls_model Optional. The base (non-bootstrapped) seminr model for full diagnostics.
#'
#' @return A cipma object containing cIPMA results, bottleneck summary, and the PLS model.
#'
#' @importFrom NCA nca_analysis
#' @importFrom stats pnorm sd
#' @export
cipma <- function(model, target_construct, data, scales, target_level = 85, nca_rep = 10000, pls_model = NULL) {

  if (is.null(model$boot_paths)) {
    stop("The seminr model must be bootstrapped to obtain CIs and p-values.")
  }

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

  nca_data <- lv_df[, c(predictors, target_construct)]

  nca_res_sig <- NCA::nca_analysis(
    data = nca_data,
    x = predictors,
    y = target_construct,
    ceilings = "ce_fdh",
    corner = 1,
    test.rep = nca_rep,
    steps = 10
  )

  nca_stats <- data.frame(
    Construct = predictors,
    NCA_d = NA,
    NCA_p_value = NA,
    Is_Necessary = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(nca_stats)) {
    pred <- nca_stats$Construct[i]
    if (!is.null(nca_res_sig$tests[[pred]])) {
      nca_stats$NCA_d[i] <- nca_res_sig$tests[[pred]]$ce_fdh$observed
      nca_stats$NCA_p_value[i] <- nca_res_sig$tests[[pred]]$ce_fdh$p_value
    } else {
      nca_stats$NCA_d[i] <- nca_res_sig$summaries[[pred]]$params[2]
    }

    d_val <- nca_stats$NCA_d[i]
    p_val <- nca_stats$NCA_p_value[i]

    if (!is.na(d_val) && !is.na(p_val)) {
      nca_stats$Is_Necessary[i] <- (d_val >= 0.10 & p_val < 0.05)
    }
  }

  bottleneck_results <- data.frame(
    Construct = predictors,
    Required_Level = NA,
    Fail_Percentage = NA,
    Fail_Count = NA,
    stringsAsFactors = FALSE
  )

  eps <- .Machine$double.eps^0.5

  for (i in 1:nrow(bottleneck_results)) {
    pred <- bottleneck_results$Construct[i]

    model_b <- NCA::nca_analysis(
      data = nca_data,
      x = pred,
      y = target_construct,
      ceilings = "ce_fdh",
      corner = 1,
      bottleneck.y = "actual",
      bottleneck.x = "actual",
      steps = c(target_level, NA)
    )

    full_bn_table <- model_b$bottlenecks$ce_fdh

    if (target_construct %in% colnames(full_bn_table)) {
      y_vals <- as.numeric(as.character(full_bn_table[[target_construct]]))
    } else {
      y_vals <- as.numeric(as.character(full_bn_table[, 1]))
    }

    target_idx <- which(abs(y_vals - target_level) < 0.001)[1]

    if (is.na(target_idx)) {
      raw_threshold <- NA
    } else {
      raw_threshold <- suppressWarnings(as.numeric(as.character(full_bn_table[target_idx, pred])))
    }

    bottleneck_results$Required_Level[i] <- raw_threshold

    if (is.na(raw_threshold)) {
      count_below <- 0
    } else {
      count_below <- sum(lv_df[[pred]] < (raw_threshold - eps), na.rm = TRUE)
    }

    bottleneck_results$Fail_Count[i] <- count_below
    bottleneck_results$Fail_Percentage[i] <- (count_below / nrow(lv_df)) * 100
  }

  final_table_A <- merge(pls_res, nca_stats, by = "Construct")
  final_table_A <- merge(final_table_A,
                         bottleneck_results[, c("Construct", "Fail_Percentage")],
                         by = "Construct")

  nca_full <- NCA::nca_analysis(
    data = nca_data,
    x = predictors,
    y = target_construct,
    ceilings = "ce_fdh",
    corner = 1,
    bottleneck.y = "actual",
    bottleneck.x = "actual",
    steps = 20
  )

  output <- list(
    cIPMA_results = final_table_A,
    Bottleneck_summary = bottleneck_results,
    Bottleneck_full = nca_full$bottlenecks$ce_fdh,
    NCA_summary = nca_res_sig,
    Target_Level = target_level,
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
