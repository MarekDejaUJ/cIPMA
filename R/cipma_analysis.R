#' Run Combined Importance-Performance Map Analysis (cIPMA)
#'
#' This function integrates PLS-SEM results with Necessary Condition Analysis (NCA).
#'
#' @param model A bootstrapped seminr model object.
#' @param target_construct String. The name of the dependent variable.
#' @param data The raw data frame used for the model.
#' @param scales A named list of theoretical min/max for indicators (e.g., list(item1=c(1,5))).
#' @param target_level Numeric. The desired level of the target construct (0-100 scale). Default is 85.
#' @param nca_rep Integer. Number of permutations for NCA. Default 10000.
#'
#' @importFrom NCA nca_analysis
#' @importFrom stats pnorm sd
#' @export
cipma <- function(model, target_construct, data, scales, target_level = 85, nca_rep = 10000) {

  # --- 1. Extract PLS-SEM Info (Importance) ---

  if (is.null(model$boot_paths)) {
    stop("The seminr model must be bootstrapped (use bootstrap_model()) to obtain CIs and p-values.")
  }

  smry <- summary(model)
  boot_total <- smry$bootstrapped_total_paths

  # Regex to find rows ending with " -> target_construct"
  pattern <- paste0("\\s*->\\s*", target_construct, "$")
  target_rows_idx <- grep(pattern, rownames(boot_total))

  if (length(target_rows_idx) == 0) {
    stop(paste0("Target construct '", target_construct, "' not found as an outcome. Check spelling."))
  }

  target_matrix <- boot_total[target_rows_idx, , drop = FALSE]
  predictors <- sub(pattern, "", rownames(target_matrix))

  pls_res <- data.frame(
    Construct = predictors,
    Importance = NA,
    Performance = NA,
    PLS_CI_Low = NA,
    PLS_CI_High = NA,
    PLS_p_value = NA,
    stringsAsFactors = FALSE
  )

  # --- 2. Calculate Rescaled Performance Scores ---

  # A) Rescale indicators to 0-100
  rescaled_data <- data
  for (item in names(scales)) {
    if(item %in% names(data)) {
      rescaled_data[[item]] <- rescale_0_100(data[[item]], scales[[item]][1], scales[[item]][2])
    }
  }

  # B) Calculate Unstandardized Weights
  outer_weights_std <- model$outer_weights
  unstandardized_weights <- outer_weights_std
  indicator_names <- rownames(outer_weights_std)

  for (indicator in indicator_names) {
    if (indicator %in% names(data)) {
      sd_indicator <- sd(data[[indicator]], na.rm = TRUE)
      if (!is.na(sd_indicator) && sd_indicator > 0) {
        unstandardized_weights[indicator, ] <- outer_weights_std[indicator, ] / sd_indicator
      }
    }
  }

  # C) Normalize Unstandardized Weights
  normalized_weights <- unstandardized_weights
  for (constr in colnames(unstandardized_weights)) {
    w_vec <- unstandardized_weights[, constr]
    total_weight <- sum(abs(w_vec), na.rm = TRUE)
    if (total_weight > 0) {
      normalized_weights[, constr] <- w_vec / total_weight
    } else {
      normalized_weights[, constr] <- 0
    }
  }

  # D) Calculate Scores
  valid_indicators <- indicator_names[indicator_names %in% colnames(rescaled_data)]
  rescaled_indicator_matrix <- as.matrix(rescaled_data[, valid_indicators])
  valid_norm_weights <- normalized_weights[valid_indicators, ]

  lv_scores_rescaled <- rescaled_indicator_matrix %*% valid_norm_weights
  lv_df <- as.data.frame(lv_scores_rescaled)

  # --- Fill PLS Results Loop ---
  for (i in 1:nrow(pls_res)) {
    pred_name <- pls_res$Construct[i]
    row_key <- rownames(target_matrix)[grep(paste0("^", pred_name, "\\s*->"), rownames(target_matrix))]

    if (length(row_key) > 0) {
      row_key <- row_key[1]
      pls_res$Importance[i] <- target_matrix[row_key, "Original Est."]

      if("2.5%" %in% colnames(target_matrix)) {
        pls_res$PLS_CI_Low[i] <- target_matrix[row_key, "2.5%"]
        pls_res$PLS_CI_High[i] <- target_matrix[row_key, "97.5%"]
      }

      cols <- colnames(target_matrix)
      if ("p" %in% cols) {
        pls_res$PLS_p_value[i] <- target_matrix[row_key, "p"]
      } else if ("T Stat." %in% cols) {
        t_val <- target_matrix[row_key, "T Stat."]
        pls_res$PLS_p_value[i] <- 2 * (1 - pnorm(abs(t_val)))
      } else {
        if(!is.na(pls_res$PLS_CI_Low[i])) {
          is_sig <- (pls_res$PLS_CI_Low[i] > 0 || pls_res$PLS_CI_High[i] < 0)
          pls_res$PLS_p_value[i] <- ifelse(is_sig, 0.049, 0.5)
        }
      }
    }

    if (pred_name %in% colnames(lv_df)) {
      pls_res$Performance[i] <- mean(lv_df[[pred_name]], na.rm = TRUE)
    }
  }

  # --- 3. Run NCA (Necessity) ---
  nca_cols <- c(predictors, target_construct)
  missing_cols <- nca_cols[!nca_cols %in% colnames(lv_df)]
  if(length(missing_cols) > 0) stop(paste("Missing constructs:", paste(missing_cols, collapse=", ")))

  nca_data <- lv_df[, nca_cols]

  nca_res <- NCA::nca_analysis(
    data = nca_data,
    x = predictors,
    y = target_construct,
    ceilings = "ce_fdh",
    corner = 1,
    test.rep = nca_rep,
    steps = 20
  )

  nca_stats <- data.frame(
    Construct = predictors,
    NCA_d = NA,
    NCA_p_value = NA,
    Is_Necessary = FALSE,
    cases_below_req = NA,
    percent_below_req = NA,
    stringsAsFactors = FALSE
  )

  # --- 4. Extract Bottlenecks (FIXED) ---
  bottleneck_full <- nca_res$bottlenecks$ce_fdh

  # FIX: Determine Y-values correctly by looking at the column, not rownames
  if (target_construct %in% colnames(bottleneck_full)) {
    y_values <- as.numeric(bottleneck_full[[target_construct]])
  } else {
    # Fallback: NCA often puts Y in the first column
    y_values <- as.numeric(bottleneck_full[, 1])
  }

  # Find the row closest to target_level
  target_row_idx <- which.min(abs(y_values - target_level))
  raw_thresholds <- bottleneck_full[target_row_idx, predictors, drop = FALSE]

  # Snapping Logic: Find closest observed value
  bottleneck_thresholds <- sapply(predictors, function(c) {
    reported_val <- suppressWarnings(as.numeric(as.character(raw_thresholds[[c]])))
    if (is.na(reported_val)) return(NA)

    vals <- lv_df[[c]]
    diffs <- abs(vals - reported_val)
    closest_val <- vals[which.min(diffs)]
    return(closest_val)
  })

  for (i in 1:nrow(nca_stats)) {
    pred <- nca_stats$Construct[i]

    if (!is.null(nca_res$tests[[pred]])) {
      nca_stats$NCA_d[i] <- nca_res$tests[[pred]]$ce_fdh$observed
      nca_stats$NCA_p_value[i] <- nca_res$tests[[pred]]$ce_fdh$p_value
    } else {
      nca_stats$NCA_d[i] <- nca_res$summaries[[pred]]$params[2]
      nca_stats$NCA_p_value[i] <- NA
    }

    d_val <- nca_stats$NCA_d[i]
    p_val <- nca_stats$NCA_p_value[i]
    if(!is.na(d_val) && !is.na(p_val)) {
      nca_stats$Is_Necessary[i] <- (d_val >= 0.10 & p_val < 0.05)
    } else {
      nca_stats$Is_Necessary[i] <- FALSE
    }

    # Bottleneck Counting using Snapped Thresholds
    thr <- bottleneck_thresholds[[pred]]

    if (is.na(thr)) {
      nca_stats$cases_below_req[i] <- 0
      nca_stats$percent_below_req[i] <- 0
    } else {
      below_count <- sum(lv_df[[pred]] < thr, na.rm = TRUE)
      nca_stats$cases_below_req[i] <- below_count
      nca_stats$percent_below_req[i] <- (below_count / nrow(lv_df)) * 100
    }
  }

  final_table_A <- merge(pls_res, nca_stats, by = "Construct")

  table_B <- data.frame(
    Construct = predictors,
    Required_Level_for_Target = bottleneck_thresholds,
    Fail_Percentage = nca_stats$percent_below_req,
    Fail_Count = nca_stats$cases_below_req
  )

  output <- list(
    cIPMA_results = final_table_A,
    Bottleneck_summary = table_B,
    Target_Level = target_level,
    Predictors = predictors,
    Dependent = target_construct,
    Data_Rescaled = lv_df
  )

  class(output) <- "cipma"
  return(output)
}
