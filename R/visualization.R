#' Plot cIPMA Results
#'
#' Visualizes the Combined Importance-Performance Map Analysis results using a bubble chart.
#'
#' @param x A cIPMA object returned by the `cipma()` function.
#' @param ... Additional arguments (not used).
#' @return A ggplot object.
#' @import ggplot2
#' @import ggrepel
#' @export
plot.cipma <- function(x, ...) {

  # Extract the main results table
  df <- x$cIPMA_results
  target <- x$Target_Level

  # --- 1. Prepare Data for Plotting ---
  df$Label <- df$Construct
  df$ColorGroup <- ifelse(df$Is_Necessary, "Necessary", "Not Necessary")

  # Size Logic: 1 for Not Necessary, % Failure for Necessary
  df$PlotSize <- ifelse(df$Is_Necessary, df$percent_below_req, 1)

  # Calculate Means for Quadrant Lines
  imp_mean <- mean(df$Importance, na.rm = TRUE)
  perf_mean <- mean(df$Performance, na.rm = TRUE)

  # --- 2. Automatic Zoom (Smart Limits) ---

  # Y-Axis (Performance) - Logic: Zoom in on data cluster
  y_range <- range(df$Performance, na.rm = TRUE)
  y_spread <- diff(y_range)
  if (y_spread == 0) y_spread <- 10
  y_pad <- y_spread * 0.15
  y_limits <- c(max(0, y_range[1] - y_pad), min(100, y_range[2] + y_pad))

  # X-Axis (Importance) - Logic: Handle Negatives + Zoom
  x_range <- range(df$Importance, na.rm = TRUE)
  x_spread <- diff(x_range)
  if (x_spread == 0) x_spread <- 0.1 # Default width if single point
  x_pad <- x_spread * 0.15

  # Calculate limits based on actual data range (supports negative values)
  # We still anchor to 0 if data is all positive to show origin context
  x_lower <- min(0, x_range[1] - x_pad)
  x_upper <- max(0, x_range[2] + x_pad)

  # If data is purely negative, ensure we see the negative side
  if (x_range[2] < 0) x_upper <- 0

  x_limits <- c(x_lower, x_upper)

  # --- 3. Build Plot ---
  p <- ggplot(df, aes(x = Importance, y = Performance)) +
    # Quadrant Lines
    geom_vline(xintercept = imp_mean, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = perf_mean, linetype = "dashed", color = "grey50") +

    # Bubbles
    geom_point(aes(size = PlotSize, fill = ColorGroup),
               shape = 21,
               color = "black",
               stroke = 0.8,
               alpha = 0.9) +

    # Labels
    geom_text_repel(aes(label = Label),
                    size = 3.5,
                    box.padding = 0.6,
                    point.padding = 0.5,
                    max.overlaps = 20) +

    # Scales
    scale_fill_manual(values = c("Necessary" = "white", "Not Necessary" = "black")) +
    scale_size_continuous(range = c(3, 12), name = "% Cases Failing\nBottleneck") +

    # Dynamic Zoom
    coord_cartesian(xlim = x_limits, ylim = y_limits) +

    # Aesthetics
    labs(
      title = paste0("cIPMA: Target Outcome Level = ", target),
      subtitle = "White bubbles = Necessary (size indicates failure rate). Black dots = Not Necessary.",
      x = "Importance (Total Effect)",
      y = "Performance (0-100)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(face = "bold", size = 14)
    )

  return(p)
}

utils::globalVariables(c("Importance", "Performance", "PlotSize", "ColorGroup", "Label", "percent_below_req", "Is_Necessary", "Construct"))
