#' Plot PCA Cumulative Explained Variance
#'
#' Visualize the cumulative explained variance of principal components (PCA)
#' and annotate the 90% and 95% thresholds to help decide the number of PCs
#' to retain for downstream analysis.
#'
#' @param obj A Seurat object. Must have the PCA reduction computed.
#' @param reduction Character. Name of the reduction to use. Default is "pca".
#' @param max_pcs Integer. Maximum number of PCs to display. Default is 50.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom tibble enframe
#' @importFrom scales percent
#' @importFrom magrittr %>%

plot_pca_variance <- function(obj, reduction = "pca", max_pcs = 50) {
  message(paste0(get_icon("step"), "[plot_pca_variance] Plotting..."))
  # 1. 检查指定的 PCA 降维是否存在
  if (!reduction %in% names(obj@reductions)) {
    stop(paste0(get_icon("error"), "[plot_pca_variance] Reduction '", reduction, "' not found in the Seurat object. Please run RunPCA first."), call. = FALSE)
  }

  # 2. 提取标准差并计算累积解释方差
  # Seurat 对象结构: obj@reductions$pca@stdev
  stdev <- obj@reductions[[reduction]]@stdev

  # 快速检查：是否有足够的 PC 数量
  if (length(stdev) < 2) {
    stop(paste0(get_icon("error"), "[plot_pca_variance] Not enough PCs found to plot."), call. = FALSE)
  }

  # 3. 逐个 PC 计算累积解释方差
  plot_data <- stdev^2 %>%
    (function(variance) cumsum(variance) / sum(variance)) %>%
    tibble::enframe(name = "PC", value = "Cumulative_Variance") %>%
    dplyr::filter(PC <= max_pcs)

  # 计算达到 90%/95% 的阈值所需最小 PC（若未达到则为 NA）
  pc_90_data <- plot_data %>% dplyr::filter(Cumulative_Variance >= 0.90)
  pc_for_90 <- if (nrow(pc_90_data) > 0) min(pc_90_data$PC) else NA

  pc_95_data <- plot_data %>% dplyr::filter(Cumulative_Variance >= 0.95)
  pc_for_95 <- if (nrow(pc_95_data) > 0) min(pc_95_data$PC) else NA

  # 4. 绘图
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

  # 若阈值存在则添加参考线与标签
  if (!is.na(pc_for_90)) {
    p <- p +
      geom_hline(yintercept = 0.90, linetype = "dashed", color = "#D55E00") +
      geom_vline(xintercept = pc_for_90, linetype = "dotted", color = "#D55E00") +
      annotate("text",
        x = pc_for_90 + 1.5, y = 0.88,
        label = paste(pc_for_90, "PCs for 90%"), color = "#D55E00", hjust = 0
      )
  }

  if (!is.na(pc_for_95)) {
    p <- p +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "#009E73") +
      geom_vline(xintercept = pc_for_95, linetype = "dotted", color = "#009E73") +
      annotate("text",
        x = pc_for_95 + 1.5, y = 0.93,
        label = paste(pc_for_95, "PCs for 95%"), color = "#009E73", hjust = 0
      )
  }

  message(paste0(get_icon("completed"), "[plot_pca_variance] Completed."))
  return(p)
}
