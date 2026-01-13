#' Rescale Data to 0-100 Range
#'
#' This helper function rescales vector data based on theoretical minimum and maximums.
#' @param x Numeric vector.
#' @param min_val Theoretical minimum of the scale.
#' @param max_val Theoretical maximum of the scale.
#' @return A numeric vector rescaled between 0 and 100.
#' @keywords internal
rescale_0_100 <- function(x, min_val, max_val) {
  if (is.na(min_val) || is.na(max_val)) return(x)
  return(((x - min_val) / (max_val - min_val)) * 100)
}

#' Get Rescaled Latent Variable Scores
#'
#' Extracts latent variable scores from a seminr model and rescales them to 0-100.
#'
#' @param model A seminr model object.
#' @param scales A named list where keys are indicator names and values are c(min, max).
#'
#' @return A data frame of rescaled latent variable scores.
#'
#' @importFrom stats sd
#' @export
get_rescaled_scores <- function(model, scales) {

  construct_scores <- model$construct_scores
  construct_names <- colnames(construct_scores)
  rescaled_scores <- as.data.frame(construct_scores)

  for (constr in construct_names) {
    mm <- model$measurement_model
    items <- mm[mm[, "construct"] == constr, "measurement"]

    if (length(items) > 0) {
      first_item <- items[1]
      if (first_item %in% names(scales)) {
        min_v <- scales[[first_item]][1]
        max_v <- scales[[first_item]][2]
        if (!is.na(min_v) && !is.na(max_v)) {
          rescaled_scores[[constr]] <- ((construct_scores[, constr] - min_v) / (max_v - min_v)) * 100
        }
      }
    }
  }

  return(rescaled_scores)
}

#' Convert bnlearn Arcs to seminr Paths
#'
#' @param learned_graph A bnlearn object (bn).
#' @return A seminr structural model object.
#' @keywords internal
arcs_to_seminr <- function(learned_graph) {
  if (!requireNamespace("bnlearn", quietly = TRUE)) {
    stop("Package 'bnlearn' is required for this function.")
  }

  arc_matrix <- bnlearn::arcs(learned_graph)

  if (nrow(arc_matrix) == 0) {
    return(NULL)
  }

  targets <- unique(arc_matrix[, "to"])
  sm_list <- list()

  for (t in targets) {
    sources <- arc_matrix[arc_matrix[, "to"] == t, "from"]
    # Create the path using seminr syntax
    sm_list[[length(sm_list) + 1]] <- seminr::paths(from = sources, to = t)
  }

  # Combine all paths into a relationship object
  return(do.call(seminr::relationships, sm_list))
}
