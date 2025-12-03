#' Run Combined Importance-Performance Map Analysis (cIPMA) - Corrected
#'
#' @param model A bootstrapped seminr model object.
#' @param target_construct String. The name of the dependent variable.
#' @param data The raw data frame used for the model.
#' @param scales A named list of theoretical min/max for indicators (e.g., list(item1=c(1,5))).
#'               Used to normalize the final LV scores.
#' @param target_level Numeric. The desired level of the target construct (0-100 scale). Default is 85.
#' @param nca_rep Integer. Number of permutations for NCA. Default 10000.
#'
#' @importFrom NCA nca_analysis
#' @importFrom stats pnorm sd
#' @export
cipma <- function(model, target_construct, data, scales, target_level = 85, nca_rep = 10000) {

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

  if(!is.null(model$sdData)) {
    std_weights_df$SD <- model$sdData[as.character(std_weights_df$Indicator)]
  } else {
    std_weights_df$SD <- apply(data[, as.character(std_weights_df$Indicator)], 2, stats::sd, na.rm=TRUE)
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

  for(i in 1:nrow(nca_stats)){
    pred <- nca_stats$Construct[i]
    if (!is.null(nca_res_sig$tests[[pred]])) {
      nca_stats$NCA_d[i] <- nca_res_sig$tests[[pred]]$ce_fdh$observed
      nca_stats$NCA_p_value[i] <- nca_res_sig$tests[[pred]]$ce_fdh$p_value
    } else {
      nca_stats$NCA_d[i] <- nca_res_sig$summaries[[pred]]$params[2]
    }

    d_val <- nca_stats$NCA_d[i]
    p_val <- nca_stats$NCA_p_value[i]

    if(!is.na(d_val) && !is.na(p_val)) {
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
      y_vals <- as.numeric(as.character(full_bn_table[,1]))
    }

    target_idx <- which(abs(y_vals - target_level) < 0.001)[1]

    if(is.na(target_idx)) {
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

  output <- list(
    cIPMA_results = final_table_A,
    Bottleneck_summary = bottleneck_results,
    Target_Level = target_level,
    Predictors = predictors,
    Dependent = target_construct,
    Data_Rescaled = lv_df
  )

  class(output) <- "cipma"
  return(output)
}
