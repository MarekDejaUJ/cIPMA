#' Run Data-Driven cIPMA with Causal Discovery
#'
#' This function performs causal discovery using constraint-based, score-based,
#' or hybrid algorithms to learn the structural model from data,
#' then runs cIPMA on the discovered structure. This exploratory approach is useful
#' when theoretical guidance is limited or when researchers want to compare
#' theory-driven and data-driven models.
#'
#' @details
#' This function implements the causal discovery approach for PLS-SEM as discussed
#' by Richter et al. (2020), who argue that combining causal discovery algorithms
#' with PLS-SEM can be valuable for:
#' \itemize{
#'   \item Exploratory research where theory is underdeveloped
#'   \item Model comparison (theory-driven vs. data-driven)
#'   \item Identifying potentially missing paths in existing models
#'   \item Cross-validation of theoretical assumptions
#' }
#'
#' The function uses proxy scores (standardized mean of indicators) for structure
#' learning, which is consistent with the composite-based nature of PLS-SEM.
#' After structure learning, proper PLS-SEM estimation is performed on the
#' original indicators.
#'
#' @section Available Algorithms:
#' \describe{
#'   \item{Score-based (fast)}{
#'     \itemize{
#'       \item \code{"hc"}: Hill-Climbing (default) - greedy search, fastest
#'       \item \code{"tabu"}: Tabu Search - avoids local optima
#'     }
#'   }
#'   \item{Constraint-based}{
#'     \itemize{
#'       \item \code{"pc"}: PC algorithm (stable version)
#'       \item \code{"gs"}: Grow-Shrink
#'       \item \code{"iamb"}: Incremental Association Markov Blanket
#'       \item \code{"fast.iamb"}: Fast IAMB
#'       \item \code{"inter.iamb"}: Interleaved IAMB
#'       \item \code{"mmpc"}: Max-Min Parents and Children
#'       \item \code{"hpc"}: Hybrid Parents and Children
#'     }
#'   }
#'   \item{Hybrid (recommended for robustness)}{
#'     \itemize{
#'       \item \code{"mmhc"}: Max-Min Hill-Climbing
#'       \item \code{"h2pc"}: Hybrid HPC
#'       \item \code{"rsmax2"}: General 2-phase Restricted Maximization
#'     }
#'   }
#' }
#'
#' @section Runtime Note:
#' The main computational cost comes from PLS-SEM bootstrapping (Step 4), not
#' structure learning. Consider reducing \code{nboot} (e.g., 500) for initial
#' exploration, then increase for final results.
#'
#' @section Caution:
#' Results from causal discovery should be interpreted as exploratory hypotheses
#' requiring theoretical justification and replication. The discovered structure
#' represents statistical associations that may reflect causal relationships but
#' cannot establish causality from observational data alone.
#'
#' @param data The raw data frame containing indicator variables.
#' @param measurement_model A seminr measurement model object (defined with constructs()).
#' @param target_construct String. The name of the dependent variable for cIPMA.
#' @param scales A named list of theoretical min/max for each indicator.
#' @param target_level Numeric. Desired target construct level (0-100). Default is 85.
#' @param nca_rep Integer. Number of permutations for NCA significance. Default 10000.
#' @param algorithm Character. Causal discovery algorithm. Options:
#'   \itemize{
#'     \item Score-based: \code{"hc"} (default), \code{"tabu"}
#'     \item Constraint-based: \code{"pc"}, \code{"gs"}, \code{"iamb"}, \code{"fast.iamb"},
#'           \code{"inter.iamb"}, \code{"mmpc"}, \code{"hpc"}
#'     \item Hybrid: \code{"mmhc"}, \code{"h2pc"}, \code{"rsmax2"}
#'   }
#' @param blacklist Optional. A data frame with columns 'from' and 'to' specifying
#'   forbidden arcs (e.g., based on temporal ordering or theory).
#' @param whitelist Optional. A data frame with columns 'from' and 'to' specifying
#'   required arcs (e.g., well-established theoretical relationships).
#' @param algorithm.args Optional. A list of additional arguments to pass to the
#'   bnlearn algorithm (e.g., \code{list(restart = 10)} for hc, \code{list(alpha = 0.01)}
#'   for constraint-based algorithms).
#' @param nboot Integer. Number of bootstrap samples for PLS-SEM. Default 1000.
#'   Use lower values (e.g., 500) for faster exploratory runs.
#'
#' @return A cipma object with additional elements:
#' \itemize{
#'   \item \code{learned_graph}: The bnlearn network object
#'   \item \code{discovered_arcs}: Data frame of discovered relationships
#'   \item \code{algorithm}: The algorithm used for discovery
#' }
#'
#' @references
#' Richter, N. F., Schubring, S., Hauff, S., Ringle, C. M., & Sarstedt, M. (2020).
#' When predictors of outcomes are necessary: Guidelines for the combined use of
#' PLS-SEM and NCA. Industrial Management & Data Systems, 120(12), 2243-2267.
#'
#' Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R Package.
#' Journal of Statistical Software, 35(3), 1-22.
#'
#' @examples
#' \dontrun{
#' # Example with whitelist (required arcs) and blacklist (forbidden arcs)
#'
#' # Whitelist: Theory strongly supports these paths
#' whitelist <- data.frame(
#'   from = c("Quality", "Satisfaction"),
#'   to = c("Satisfaction", "Loyalty")
#' )
#'
#' # Blacklist: These paths are theoretically impossible
#' blacklist <- data.frame(
#'   from = c("Loyalty", "Loyalty"),
#'   to = c("Quality", "Satisfaction")
#' )
#'
#' # Run with hybrid algorithm (recommended)
#' result <- discovery_cipma(
#'   data = survey_data,
#'   measurement_model = mm,
#'   target_construct = "Loyalty",
#'   scales = scales,
#'   algorithm = "mmhc",
#'   whitelist = whitelist,
#'   blacklist = blacklist,
#'   nboot = 500  # Lower for faster exploration
#' )
#'
#' # View discovered structure
#' print(result$discovered_arcs)
#' }
#'
#' @importFrom stats sd
#' @export
discovery_cipma <- function(data, measurement_model, target_construct, scales,
                            target_level = 85, nca_rep = 10000,
                            algorithm = "hc", blacklist = NULL, whitelist = NULL,
                            algorithm.args = list(), nboot = 1000) {

  if (!requireNamespace("bnlearn", quietly = TRUE)) {
    stop("Package 'bnlearn' is required for discovery_cipma(). ",
         "Please install it with: install.packages('bnlearn')")
  }

  score_based <- c("hc", "tabu")
  constraint_based <- c("pc", "gs", "iamb", "fast.iamb", "inter.iamb", "mmpc", "hpc")
  hybrid <- c("mmhc", "h2pc", "rsmax2")

  all_algorithms <- c(score_based, constraint_based, hybrid)
  if (!algorithm %in% all_algorithms) {
    stop("Unknown algorithm '", algorithm, "'. Choose from: ",
         paste(all_algorithms, collapse = ", "))
  }

  message("=== Data-Driven cIPMA (Causal Discovery) ===\n")
  message("--- Step 1: Calculating Construct Scores for Structure Learning ---")

  construct_names <- sapply(measurement_model, function(x) x[[1]])

  if (length(construct_names) >= 2) {
    path_list <- list()
    for (i in 2:length(construct_names)) {
      path_list[[i-1]] <- seminr::paths(from = construct_names[1], to = construct_names[i])
    }
    initial_sm <- do.call(seminr::relationships, path_list)
  } else {
    stop("At least 2 constructs are required for structure learning.")
  }

  temp_pls <- tryCatch({
    suppressMessages(seminr::estimate_pls(data, measurement_model, initial_sm))
  }, error = function(e) NULL)

  if (!is.null(temp_pls) && ncol(temp_pls$construct_scores) == length(construct_names)) {
    scores_df <- as.data.frame(temp_pls$construct_scores)
    message("    Using PLS-SEM construct scores (recommended)")
  } else {
    scores_df <- data.frame(matrix(ncol = length(construct_names), nrow = nrow(data)))
    colnames(scores_df) <- construct_names

    for (i in seq_along(measurement_model)) {
      c_name  <- measurement_model[[i]][1]
      c_items <- measurement_model[[i]][-1]

      valid_items <- c_items[c_items %in% colnames(data)]
      if (length(valid_items) == 0) {
        stop(paste("No valid data columns found for construct:", c_name))
      }

      item_data <- data[, valid_items, drop = FALSE]
      item_data <- as.data.frame(lapply(item_data, as.numeric))
      scores_df[[c_name]] <- scale(rowMeans(item_data, na.rm = TRUE))
    }
    message("    Using proxy scores (mean of indicators)")
  }

  node_names <- colnames(scores_df)
  message("    Constructs: ", paste(node_names, collapse = ", "))

  if (!is.null(blacklist)) {
    bl_nodes <- unique(c(blacklist$from, blacklist$to))
    invalid_bl <- bl_nodes[!bl_nodes %in% node_names]
    if (length(invalid_bl) > 0) {
      stop("Blacklist contains unknown construct names: ", paste(invalid_bl, collapse = ", "),
           "\n  Available constructs: ", paste(node_names, collapse = ", "))
    }
  }

  if (!is.null(whitelist)) {
    wl_nodes <- unique(c(whitelist$from, whitelist$to))
    invalid_wl <- wl_nodes[!wl_nodes %in% node_names]
    if (length(invalid_wl) > 0) {
      stop("Whitelist contains unknown construct names: ", paste(invalid_wl, collapse = ", "),
           "\n  Available constructs: ", paste(node_names, collapse = ", "))
    }
  }

  algo_type <- ifelse(algorithm %in% score_based, "Score-based",
                      ifelse(algorithm %in% constraint_based, "Constraint-based", "Hybrid"))

  message("--- Step 2: Learning Causal Structure ---")
  message("    Algorithm: ", algorithm, " (", algo_type, ")")

  base_args <- list(x = scores_df)
  if (!is.null(blacklist)) base_args$blacklist <- blacklist
  if (!is.null(whitelist)) base_args$whitelist <- whitelist

  all_args <- c(base_args, algorithm.args)

  learned_graph <- switch(algorithm,
                          "hc" = do.call(bnlearn::hc, all_args),
                          "tabu" = do.call(bnlearn::tabu, all_args),
                          "pc" = do.call(bnlearn::pc.stable, all_args),
                          "gs" = do.call(bnlearn::gs, all_args),
                          "iamb" = do.call(bnlearn::iamb, all_args),
                          "fast.iamb" = do.call(bnlearn::fast.iamb, all_args),
                          "inter.iamb" = do.call(bnlearn::inter.iamb, all_args),
                          "mmpc" = do.call(bnlearn::mmpc, all_args),
                          "hpc" = do.call(bnlearn::hpc, all_args),
                          "mmhc" = do.call(bnlearn::mmhc, all_args),
                          "h2pc" = do.call(bnlearn::h2pc, all_args),
                          "rsmax2" = do.call(bnlearn::rsmax2, all_args)
  )

  arc_matrix <- bnlearn::arcs(learned_graph)

  if (nrow(arc_matrix) == 0) {
    stop("Causal discovery found no relationships. Consider:\n",
         "  - Adjusting algorithm parameters (e.g., alpha for constraint-based)\n",
         "  - Using a different algorithm\n",
         "  - Checking data quality and sample size")
  }

  message("    Discovered ", nrow(arc_matrix), " arcs:")
  for (i in seq_len(nrow(arc_matrix))) {
    message("      ", arc_matrix[i, "from"], " -> ", arc_matrix[i, "to"])
  }

  if (!is.null(whitelist)) {
    message("    Whitelist (required): ", nrow(whitelist), " arcs")
  }
  if (!is.null(blacklist)) {
    message("    Blacklist (forbidden): ", nrow(blacklist), " arcs")
  }

  targets <- unique(arc_matrix[, "to"])
  sm_list <- list()
  for (t in targets) {
    sources <- arc_matrix[arc_matrix[, "to"] == t, "from"]
    sm_list[[length(sm_list) + 1]] <- seminr::paths(from = sources, to = t)
  }

  discovered_sm <- do.call(seminr::relationships, sm_list)

  message("\n--- Step 3: Running PLS-SEM on Discovered Structure ---")

  pls_model <- seminr::estimate_pls(data, measurement_model, discovered_sm)

  message("--- Step 4: Bootstrapping Model (", nboot, " samples) ---")
  message("    (This may take a while...)")

  boot_model <- seminr::bootstrap_model(pls_model, nboot = nboot)

  message("--- Step 5: Running cIPMA Analysis ---")

  result <- cipma(
    model = boot_model,
    target_construct = target_construct,
    data = data,
    scales = scales,
    target_level = target_level,
    nca_rep = nca_rep,
    pls_model = pls_model
  )

  result$learned_graph <- learned_graph
  result$discovered_arcs <- as.data.frame(arc_matrix)
  result$algorithm <- algorithm
  result$algorithm_type <- algo_type

  message("\n=== Discovery cIPMA Complete ===")
  message("Algorithm: ", algorithm, " (", algo_type, ")")
  message("Use print(result$discovered_arcs) to view learned structure.")
  message("Use plot(result$learned_graph) to visualize the DAG.")
  message("Use check_assumptions(result) to validate model quality.")

  result
}
