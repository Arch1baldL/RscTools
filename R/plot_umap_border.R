#' Plot UMAP with Density Borders
#'
#' This function visualizes single-cell UMAP embeddings with scatter points and
#' draws density borders around specified cell groups. It automatically handles
#' UMAP axes with arrows and expands the plotting area to prevent border clipping.
#' It auto-detects coordinate names based on the provided reduction.
#'
#' @param obj A Seurat object containing UMAP coordinates.
#' @param group_points String. The metadata column name to use for coloring scatter points (e.g., "cell_cluster").
#' @param group_borders String. The metadata column name to use for drawing borders (e.g., "cell_sub_type").
#'   If NULL, defaults to `group_points`.
#' @param reduction String. The dimensional reduction to use (default is "umap").
#' @param min_cells Integer. Minimum number of cells required in a group to draw a border. Groups with fewer cells will only show scatter points. Default is 10.
#' @param expand_ratio Numeric. The ratio to expand the bounding box for density calculation to avoid clipping. Default is 0.3 (30%).
#' @param pt_size Numeric. Size of the scatter points. Default is 0.5.
#' @param line_width Numeric. Width of the border lines. Default is 1.2.
#' @param label_size Numeric. Font size for UMAP axis labels. Default is 5.
#' @param border_breaks Numeric. Density threshold for drawing borders (passed to stat_density_2d breaks). Default is 0.15.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData
#' @importFrom ggplot2 ggplot aes geom_point stat_density_2d annotate theme_void theme coord_fixed scale_color_discrete margin arrow unit
#' @importFrom dplyr filter group_by group_modify ungroup bind_rows %>%
#' @importFrom grid arrow unit
#' @importFrom rlang .data

plot_umap_border <- function(
  obj,
  group_points,
  group_borders = NULL,
  reduction = "umap",
  min_cells = 10,
  expand_ratio = 0.3,
  pt_size = 0.5,
  line_width = 1.2,
  label_size = 5,
  border_breaks = 0.1
) {
  message(paste0(get_icon("step"), "[plot_umap_border] Plotting..."))
  if (is.null(group_borders)) {
    group_borders <- group_points
  }
  # 1. 校验与坐标提取
  # 检查 reduction 是否存在
  if (!reduction %in% names(obj@reductions)) {
    stop(paste0(get_icon("error"), "[plot_umap_border] Reduction '", reduction, "' not found in Seurat object. Available reductions: ", paste(names(obj@reductions), collapse = ", ")), call. = FALSE)
  }

  # 动态获取坐标列名（如 "umap_1"、"UMAP_1"、"PC_1"）
  # 使用 [[ ]] 访问 DimReduc 对象并读取列名
  coord_names <- colnames(obj[[reduction]])[1:2]

  # 使用正确列名提取坐标与元数据
  umap_data <- Seurat::FetchData(obj, vars = c(coord_names, group_points, group_borders))

  # 统一内部列名（标准化为 UMAP_1/UMAP_2 以复用绘图逻辑）
  colnames(umap_data) <- c("UMAP_1", "UMAP_2", "group_for_points", "group_for_hulls")

  # 2. 过滤：仅保留满足最小细胞数的分组用于绘制边界
  group_counts <- table(umap_data$group_for_hulls)
  valid_groups <- names(group_counts[group_counts >= min_cells])

  if (length(valid_groups) == 0) {
    stop(paste0(get_icon("error"), "[plot_umap_border] All groups have fewer cells than 'min_cells'. Cannot draw borders."), call. = FALSE)
  }

  umap_data_for_contour <- umap_data %>%
    dplyr::filter(.data$group_for_hulls %in% valid_groups)

  # 3. 扩展边界（幽灵点机制）
  umap_data_expanded <- umap_data_for_contour %>%
    dplyr::group_by(.data$group_for_hulls) %>%
    dplyr::group_modify(~ {
      data <- .x
      range_x <- range(data$UMAP_1)
      range_y <- range(data$UMAP_2)
      span_x <- diff(range_x) * expand_ratio
      span_y <- diff(range_y) * expand_ratio

      ghost_points <- data.frame(
        UMAP_1 = c(range_x[1] - span_x, range_x[2] + span_x, range_x[1] - span_x, range_x[2] + span_x),
        UMAP_2 = c(range_y[1] - span_y, range_y[2] + span_y, range_y[2] + span_y, range_y[1] - span_y),
        group_for_points = NA
      )
      # dplyr::bind_rows 会利用分组上下文自动补齐缺失的分组列
      dplyr::bind_rows(data, ghost_points)
    }) %>%
    dplyr::ungroup()

  # 4. 计算坐标轴位置
  min_x <- min(umap_data$UMAP_1)
  min_y <- min(umap_data$UMAP_2)
  max_x <- max(umap_data$UMAP_1)

  arrow_len <- (max_x - min_x) * 0.15
  offset <- 1
  origin_x <- min_x - offset
  origin_y <- min_y - offset
  label_gap <- (max_x - min_x) * 0.02

  # 5. 绘图
  p <- ggplot2::ggplot(umap_data, ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)) +

    # A. 散点
    ggplot2::geom_point(ggplot2::aes(color = .data$group_for_points), size = pt_size, alpha = 0.5) +

    # B. 边界
    ggplot2::stat_density_2d(
      data = umap_data_expanded,
      ggplot2::aes(fill = .data$group_for_hulls, color = .data$group_for_hulls),
      geom = "polygon",
      contour_var = "ndensity",
      breaks = border_breaks,
      size = line_width,
      alpha = 0.2,
      show.legend = FALSE
    ) +

    # C. 箭头坐标轴
    # X 轴
    ggplot2::annotate("segment",
      x = origin_x, xend = origin_x + arrow_len,
      y = origin_y, yend = origin_y,
      arrow = ggplot2::arrow(type = "open", length = ggplot2::unit(0.1, "inches")),
      size = 1, color = "black"
    ) +
    ggplot2::annotate("text",
      x = origin_x + (arrow_len * 0.05), y = origin_y - label_gap,
      label = "UMAP_1", hjust = 0, vjust = 1,
      size = label_size, fontface = "bold"
    ) +

    # Y 轴
    ggplot2::annotate("segment",
      x = origin_x, xend = origin_x,
      y = origin_y, yend = origin_y + arrow_len,
      arrow = ggplot2::arrow(type = "open", length = ggplot2::unit(0.1, "inches")),
      size = 1, color = "black"
    ) +
    ggplot2::annotate("text",
      x = origin_x - label_gap, y = origin_y + (arrow_len * 0.05),
      label = "UMAP_2", hjust = 0, vjust = 0, angle = 90,
      size = label_size, fontface = "bold"
    ) +

    # D. 版式
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20)) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::scale_color_discrete(name = "Group")

  message(paste0(get_icon("completed"), "[plot_umap_border] Completed."))
  return(p)
}
