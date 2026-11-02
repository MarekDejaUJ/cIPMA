#' Plot cIPMA Results
#'
#' Draws the combined map, the tolerance path of the effect sizes or the
#' tolerance path of the required levels, in grey scale. On the map, white
#' bubbles mark the predictors necessary at \code{pi = 1}, with a size given by
#' the share of cases below the required level at the chosen tolerance, and
#' black dots mark the other predictors. \code{scale = "zoom"} fits the axes to
#' the results; \code{scale = "full"} shows performance on the whole 0-100 scale
#' and importance from 0 to the next tenth above the largest total effect.
#'
#' @param x A cipma object returned by \code{\link{cipma}}.
#' @param type One of \code{"map"}, \code{"path"} (effect size across
#'   tolerances) or \code{"bottleneck"} (required level at the target level
#'   across tolerances).
#' @param scale One of \code{"zoom"} or \code{"full"}; used by the map.
#' @param ... Ignored.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @import ggrepel
#' @export
plot.cipma <- function(x, type = c("map", "path", "bottleneck"), scale = c("zoom", "full"), ...) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  if (type == "map") return(cipma_map(x, scale))
  cipma_path_plot(x, type)
}

cipma_map <- function(x, scale) {
  df <- x$cIPMA_results
  target <- x$Target_Level

  df$Label <- df$Construct
  necessary <- !is.na(df$Is_Necessary) & df$Is_Necessary
  df$ColorGroup <- ifelse(necessary, "Necessary", "Not Necessary")
  df$PlotSize <- ifelse(necessary, df$Fail_Percentage, 1)

  imp_mean <- mean(df$Importance, na.rm = TRUE)
  perf_mean <- mean(df$Performance, na.rm = TRUE)

  x_range <- range(df$Importance, na.rm = TRUE)
  if (scale == "full") {
    y_limits <- c(0, 100)
    x_limits <- c(min(0, floor(x_range[1] * 10) / 10), max(0.1, ceiling(x_range[2] * 10) / 10))
  } else {
    y_range <- range(df$Performance, na.rm = TRUE)
    y_spread <- diff(y_range)
    if (y_spread == 0) y_spread <- 10
    y_pad <- y_spread * 0.15
    y_limits <- c(max(0, y_range[1] - y_pad), min(100, y_range[2] + y_pad))
    x_spread <- diff(x_range)
    if (x_spread == 0) x_spread <- 0.1
    x_pad <- x_spread * 0.15
    x_lower <- min(0, x_range[1] - x_pad)
    x_upper <- max(0, x_range[2] + x_pad)
    if (x_range[2] < 0) x_upper <- 0
    x_limits <- c(x_lower, x_upper)
  }

  ggplot(df, aes(x = Importance, y = Performance)) +
    geom_vline(xintercept = imp_mean, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = perf_mean, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = PlotSize, fill = ColorGroup),
               shape = 21, color = "black", stroke = 0.8, alpha = 0.9) +
    geom_text_repel(aes(label = Label), size = 3.5, box.padding = 0.6,
                    point.padding = 0.5, max.overlaps = 20, colour = "grey10") +
    scale_fill_manual(values = c("Necessary" = "white", "Not Necessary" = "black"), name = NULL) +
    scale_size_continuous(range = c(3, 12),
                          name = sprintf("%% cases failing\nbottleneck (pi = %s)", pi_label(x$Tolerance))) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    labs(
      title = sprintf("cIPMA: %s at level %s", x$Dependent, format(target)),
      subtitle = "White bubbles: necessary at pi = 1 (size: share of cases failing the bottleneck). Black dots: not necessary.",
      x = "Importance (total effect)",
      y = "Performance (0-100)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9)
    )
}

cipma_path_plot <- function(x, type) {
  df <- x$Tolerance_Path
  shapes <- rep(c(21, 22, 24, 23, 25, 1, 0, 2, 5, 6), length.out = length(unique(df$Construct)))
  y_col <- if (type == "path") "d" else "Required_Level"
  y_lab <- if (type == "path") "Effect size d" else
    sprintf("Required level for %s >= %s", x$Dependent, format(x$Target_Level))
  title <- if (type == "path") "Tolerance path of the effect size" else "Tolerance path of the bottleneck"
  df$y <- df[[y_col]]
  pis <- sort(unique(df$pi))
  ggplot(df, aes(x = pi, y = .data$y, group = Construct, linetype = Construct, shape = Construct)) +
    geom_line(colour = "grey25", na.rm = TRUE) +
    geom_point(colour = "grey10", fill = "white", size = 2, na.rm = TRUE) +
    scale_x_reverse(breaks = pis, labels = pi_label(pis)) +
    scale_shape_manual(values = shapes) +
    labs(title = sprintf("%s: %s", title, x$Dependent), x = "Tolerance pi", y = y_lab) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), axis.line = element_line(color = "black"),
          plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(size = 7))
}
