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

  # Label: Remove "(Necessary)" text to reduce clutter; color/fill already tells the story
  df$Label <- df$Construct

  # Color Groups: Match Paper (White = Necessary, Black = Not Necessary)
  df$ColorGroup <- ifelse(df$Is_Necessary, "Necessary", "Not Necessary")

  # Size Logic:
  # If Necessary: Size = % Failing (Large bubbles = Urgency)
  # If Not Necessary: Size = Fixed small value (Standard dot)
  # This mimics Fig 4 in Hauff et al. (2024)
  df$PlotSize <- ifelse(df$Is_Necessary, df$percent_below_req, 1)

  # Calculate Means for Quadrants
  imp_mean <- mean(df$Importance, na.rm = TRUE)
  perf_mean <- mean(df$Performance, na.rm = TRUE)

  # --- 2. Automatic Zoom (Smart Limits) ---

  # Y-Axis (Performance)
  y_range <- range(df$Performance, na.rm = TRUE)
  spread <- diff(y_range)

  # SAFETY FIX: If data has no spread (single point or flat line), force a window
  if (spread == 0) spread <- 10

  padding <- spread * 0.15
  y_limits <- c(max(0, y_range[1] - padding), min(100, y_range[2] + padding))

  # X-Axis (Importance)
  x_max <- max(df$Importance, na.rm = TRUE)
  # Ensure we show 0 for context, but don't squish too much
  x_limits <- c(0, ifelse(x_max == 0, 1, x_max * 1.2))

  # --- 3. Build Plot ---
  p <- ggplot(df, aes(x = Importance, y = Performance)) +
    # Quadrant Lines
    geom_vline(xintercept = imp_mean, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = perf_mean, linetype = "dashed", color = "grey50") +

    # Bubbles
    # Note: We map size to 'PlotSize' calculated above
    geom_point(aes(size = PlotSize, fill = ColorGroup),
               shape = 21,       # Allows fill (white/black) and color (black outline)
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
    # Match Paper: White bubbles for necessary, Black dots for not necessary
    scale_fill_manual(values = c("Necessary" = "white", "Not Necessary" = "black")) +

    # Size Scale: We set the range so '1' (Not necessary) is small (3), others are larger (up to 12)
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
