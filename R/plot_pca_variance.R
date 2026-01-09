#' Plot PCA Cumulative Explained Variance
#'
#' This function visualizes the cumulative explained variance of Principal Components (PCA).
#' It helps in determining the optimal number of PCs to retain for downstream analysis
#' by marking the 90% and 95% variance thresholds.
#'
#' @param obj A Seurat object. Must have PCA reduction computed.
#' @param reduction Character. Name of the reduction to use. Default is "pca".
#' @param max_pcs Integer. Maximum number of PCs to plot. Default is 50.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom tibble enframe
#' @importFrom scales percent
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' plot_pca_variance(seurat_obj)
#' }

plot_pca_variance <- function(obj, reduction = "pca", max_pcs = 50) {

  # 1. Check whether the specified PCA reduction exists
  if (!reduction %in% names(obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in the Seurat object. Please run RunPCA first."))
  }

  # 2. Extract standard deviations and compute cumulative explained variance
  # Seurat object structure: obj@reductions$pca@stdev
  stdev <- obj@reductions[[reduction]]@stdev

  # Quick sanity check for a sufficient number of PCs
  if (length(stdev) < 2) {
    stop("Not enough PCs found to plot.")
  }

  # 3. Compute cumulative explained variance per PC
  plot_data <- stdev^2 %>%
    (function(variance) cumsum(variance) / sum(variance)) %>%
    tibble::enframe(name = "PC", value = "Cumulative_Variance") %>%
    dplyr::filter(PC <= max_pcs)

  # Compute threshold PCs (robust in case 90%/95% are not reached)
  pc_90_data <- plot_data %>% dplyr::filter(Cumulative_Variance >= 0.90)
  pc_for_90 <- if (nrow(pc_90_data) > 0) min(pc_90_data$PC) else NA

  pc_95_data <- plot_data %>% dplyr::filter(Cumulative_Variance >= 0.95)
  pc_for_95 <- if (nrow(pc_95_data) > 0) min(pc_95_data$PC) else NA

  # 4. Visualization
  p <- ggplot(plot_data, aes(x = PC, y = Cumulative_Variance)) +
    geom_line(color = "#0072B2", linewidth = 1) +
    geom_point(color = "#0072B2", size = 2.5, alpha = 0.8) +
    labs(
      title = "PCA Cumulative Explained Variance",
      x = "Number of Principal Components",
      y = "Cumulative Explained Variance Ratio"
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(NA, 1.05), n.breaks = 6) +
    scale_x_continuous(breaks = seq(0, max_pcs, by = 5)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  # Add dynamic reference lines and labels if thresholds exist
  if (!is.na(pc_for_90)) {
    p <- p +
      geom_hline(yintercept = 0.90, linetype = "dashed", color = "#D55E00") +
      geom_vline(xintercept = pc_for_90, linetype = "dotted", color = "#D55E00") +
      annotate("text", x = pc_for_90 + 1.5, y = 0.88,
               label = paste(pc_for_90, "PCs for 90%"), color = "#D55E00", hjust = 0)
  }

  if (!is.na(pc_for_95)) {
    p <- p +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "#009E73") +
      geom_vline(xintercept = pc_for_95, linetype = "dotted", color = "#009E73") +
      annotate("text", x = pc_for_95 + 1.5, y = 0.93,
               label = paste(pc_for_95, "PCs for 95%"), color = "#009E73", hjust = 0)
  }

  return(p)
}