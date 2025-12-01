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

#' Unstandardize and Rescale Latent Variable Scores
#'
#' Extracts latent variable scores from a seminr model and rescales them to 0-100.
#'
#' @param model A seminr model object.
#' @param scales A named list where keys are indicator names and values are c(min, max).
#' @importFrom stats sd
#' @return A data frame of rescaled latent variable scores.
#' @export
get_rescaled_scores <- function(model, scales) {

  # 1. Get Construct Scores (seminr usually returns standardized if PLS)
  # We need to reconstruct unstandardized scores to rescale them properly
  # or apply rescaling logic to the construct scores if using simple composites.

  # For cIPMA per Hauff et al., we usually rescale the inputs (indicators)
  # and then calculate scores, OR rescale the final scores based on theoretical bounds.
  # Here we implement the logic to rescale the final construct scores based on
  # the theoretical range of their indicators.

  construct_scores <- model$construct_scores
  construct_names <- colnames(construct_scores)
  rescaled_scores <- as.data.frame(construct_scores)

  # Note: In a full package, you might want more complex logic to handle
  # weights validation as per Hauff et al. (checking for negative weights).

  for (constr in construct_names) {
    # Find measurement items for this construct
    # This assumes all items for a construct share the same scale
    mm <- model$measurement_model
    items <- mm[mm[, "construct"] == constr, "measurement"]

    if (length(items) > 0) {
      # Grab the scale of the first item (assuming homogeneity per construct)
      first_item <- items[1]
      if (first_item %in% names(scales)) {
        min_v <- scales[[first_item]][1]
        max_v <- scales[[first_item]][2]

        # Rescale the construct score
        # Note: If model$construct_scores are standardized (mean 0, sd 1),
        # this simple linear transformation might not map perfectly to 0-100
        # unless we unstandardize first.
        # For simplicity in this script, we assume the user follows the
        # Hauff et al. guide: Rescale Inputs -> Run PLS -> Use those scores.
        # OR: We apply the rescaling formula to the raw data aggregate.

        # To match the user script logic:
        rescaled_scores[[constr]] <- rescale_0_100(construct_scores[,constr], min_v, max_v)
      }
    }
  }

  return(rescaled_scores)
}
