#' Run Data-Driven cIPMA (Causal Discovery)
#'
#' This function performs a Causal Discovery analysis using the Hill-Climbing
#' algorithm (bnlearn) to learn the structural model from data, and then runs
#' cIPMA on the discovered structure.
#'
#' @param data The raw data frame.
#' @param measurement_model A seminr measurement model object (defined with constructs()).
#' @param target_construct String. The name of the dependent variable.
#' @param scales A named list of theoretical min/max for indicators.
#' @param target_level Numeric. Desired target construct level (0–100). Default is 85.
#' @param nca_rep Integer. Number of permutations for NCA. Default 10000.
#' @param blacklist Optional. A data frame with columns 'from' and 'to' defining
#'   forbidden arcs for bnlearn::hc.
#'
#' @return A cIPMA object (same class as cipma()), with an additional
#'   `learned_graph` element containing the bnlearn network.
#'
#' @importFrom seminr estimate_pls bootstrap_model relationships paths
#' @importFrom bnlearn hc arcs
#' @importFrom stats sd
#' @export
discovery_cipma <- function(data, measurement_model, target_construct, scales,
                            target_level = 85, nca_rep = 10000, blacklist = NULL) {

  message("--- Step 1: Calculating Proxy Scores for Structure Learning ---")

  construct_names <- sapply(measurement_model, function(x) x[[1]])
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

  message("--- Step 2: Learning Causal Structure (Hill-Climbing) ---")

  learned_graph <- bnlearn::hc(scores_df, blacklist = blacklist)
  arc_matrix    <- bnlearn::arcs(learned_graph)

  if (nrow(arc_matrix) == 0) {
    stop("Causal discovery found no relationships (no arcs). Cannot proceed with PLS.")
  }

  targets <- unique(arc_matrix[, "to"])
  sm_list <- list()
  for (t in targets) {
    sources <- arc_matrix[arc_matrix[, "to"] == t, "from"]
    sm_list[[length(sm_list) + 1]] <- seminr::paths(from = sources, to = t)
  }

  discovered_sm <- do.call(seminr::relationships, sm_list)

  message("--- Step 3: Running PLS-SEM on Discovered Structure ---")

  final_pls  <- seminr::estimate_pls(data, measurement_model, discovered_sm)

  message("--- Step 4: Bootstrapping Model ---")
  boot_final <- seminr::bootstrap_model(final_pls, nboot = 1000)

  message("--- Step 5: Running cIPMA Analysis ---")

  result <- cipma(
    model          = boot_final,
    target_construct = target_construct,
    data           = data,
    scales         = scales,
    target_level   = target_level,
    nca_rep        = nca_rep
  )

  result$learned_graph <- learned_graph

  result
}
