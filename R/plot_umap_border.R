#' Plot UMAP with Density Borders and Labels
#'
#' Visualize UMAP embeddings with scatter points, density borders, and optional
#' group labels. Handles axis arrows and expands the plotting area to prevent
#' border clipping.
#'
#' @param obj A Seurat object containing embeddings.
#' @param group_points Metadata column for coloring points.
#' @param group_borders Metadata column for borders/labels; defaults to `group_points` when NULL or missing.
#' @param reduction Dimensional reduction to use (default "umap").
#' @param min_cells Minimum cells per group to draw a border. Default 10.
#' @param expand_ratio Expand ratio to avoid clipping. Default 0.3.
#' @param pt_size Point size. Default 0.5.
#' @param line_width Border line width. Default 1.2.
#' @param label_size Axis label size. Default 5.
#' @param show_labels Logical; show group labels. Default TRUE.
#' @param text_size Label text size. Default 4.
#' @param smoothness Grid density for stat_density_2d; larger gives smoother contours. Default 300.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData Embeddings
#' @importFrom ggplot2 ggplot aes geom_point stat_density_2d annotate theme_void theme coord_fixed scale_color_discrete margin arrow unit
#' @importFrom dplyr filter group_by group_modify ungroup bind_rows summarise %>%
#' @importFrom grid arrow unit
#' @importFrom rlang .data
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom stats median

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
  show_labels = TRUE,
  text_size = 4,
  smoothness = 300
) {
  message(paste0(get_icon("step"), "[plot_umap_border] Plotting..."))
  if (is.null(group_borders)) {
    group_borders <- group_points
  }
  # 参数已取消：以下为内部默认设置（如需调整可改此处）
  border_breaks <- 0.1 # 原参数 border_breaks 的默认值

  # 1. 校验与坐标提取
  if (!reduction %in% names(obj@reductions)) {
    stop(paste0(get_icon("error"), "[plot_umap_border] Reduction '", reduction, "' not found in Seurat object. Available reductions: ", paste(names(obj@reductions), collapse = ", ")), call. = FALSE)
  }

  coord_names <- colnames(obj[[reduction]])[1:2]

  # 元数据列校验
  meta_cols <- colnames(obj[[]])
  if (!group_points %in% meta_cols) {
    stop(paste0(get_icon("error"), "[plot_umap_border] Metadata column '", group_points, "' not found."), call. = FALSE)
  }
  if (!group_borders %in% meta_cols) {
    warning(paste0(get_icon("warning"), "[plot_umap_border] Metadata column '", group_borders, "' not found. Using '", group_points, "' for borders."), call. = FALSE)
    group_borders <- group_points
  }

  vars <- unique(c(coord_names, group_points, group_borders))
  df <- Seurat::FetchData(obj, vars = vars)
  umap_data <- data.frame(
    UMAP_1 = df[[coord_names[1]]],
    UMAP_2 = df[[coord_names[2]]],
    group_for_points = df[[group_points]],
    group_for_hulls = df[[group_borders]],
    check.names = FALSE
  )
  if (is.factor(umap_data$group_for_hulls)) umap_data$group_for_hulls <- droplevels(umap_data$group_for_hulls)
  if (is.factor(umap_data$group_for_points)) umap_data$group_for_points <- droplevels(umap_data$group_for_points)
  umap_data <- umap_data[!is.na(umap_data$group_for_hulls), , drop = FALSE]

  # 2. 过滤：仅保留满足最小细胞数的分组用于绘制边界
  group_counts <- table(umap_data$group_for_hulls)
  valid_groups <- names(group_counts[group_counts >= min_cells])

  if (length(valid_groups) == 0) {
    stop(paste0(get_icon("error"), "[plot_umap_border] All groups have fewer cells than 'min_cells'. Cannot draw borders."), call. = FALSE)
  }

  umap_data_for_contour <- umap_data %>%
    dplyr::filter(!is.na(.data$group_for_hulls)) %>%
    dplyr::filter(.data$group_for_hulls %in% valid_groups)

  # 3. 扩展边界（幽灵点机制）
  umap_data_expanded <- umap_data_for_contour %>%
    dplyr::group_split(.data$group_for_hulls) %>%
    lapply(function(data) {
      if (nrow(data) == 0) {
        return(NULL)
      }
      range_x <- range(data$UMAP_1)
      range_y <- range(data$UMAP_2)
      span_x <- diff(range_x) * expand_ratio
      span_y <- diff(range_y) * expand_ratio
      grp_val <- data$group_for_hulls[1]
      if (length(grp_val) == 0) grp_val <- NA_character_

      ghost_points <- data.frame(
        UMAP_1 = c(range_x[1] - span_x, range_x[2] + span_x, range_x[1] - span_x, range_x[2] + span_x),
        UMAP_2 = c(range_y[1] - span_y, range_y[2] + span_y, range_y[2] + span_y, range_y[1] - span_y),
        group_for_points = rep(NA_character_, 4),
        group_for_hulls = rep(grp_val, 4)
      )
      dplyr::bind_rows(data, ghost_points)
    }) %>%
    dplyr::bind_rows()

  # 4. 计算坐标轴位置
  min_x <- min(umap_data$UMAP_1)
  min_y <- min(umap_data$UMAP_2)
  max_x <- max(umap_data$UMAP_1)

  arrow_len <- (max_x - min_x) * 0.15
  offset <- 1
  origin_x <- min_x - offset
  origin_y <- min_y - offset
  label_gap <- (max_x - min_x) * 0.02

  # 5. 标签位置
  if (isTRUE(show_labels)) {
    label_data <- umap_data_for_contour %>%
      dplyr::group_by(group_for_hulls) %>%
      dplyr::summarise(
        UMAP_1 = stats::median(UMAP_1),
        UMAP_2 = stats::median(UMAP_2),
        .groups = "drop"
      )
  }

  # 6. 绘图
  p <- ggplot2::ggplot(umap_data, ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)) +

    # A. 散点
    ggplot2::geom_point(ggplot2::aes(color = .data$group_for_points), size = pt_size, alpha = 0.5) +

    # B. 边界主线层（虚线、平滑）
    ggplot2::stat_density_2d(
      data = umap_data_expanded,
      ggplot2::aes(color = .data$group_for_hulls),
      geom = "contour",
      contour_var = "ndensity",
      breaks = border_breaks, # 使用内部默认阈值
      size = line_width,
      linetype = "dashed",
      n = smoothness,
      show.legend = FALSE
    )

  if (isTRUE(show_labels)) {
    p <- p + ggrepel::geom_label_repel(
      data = label_data,
      ggplot2::aes(x = UMAP_1, y = UMAP_2, label = group_for_hulls, color = group_for_hulls),
      fill = "white",
      label.size = 0.6,
      label.r = ggplot2::unit(0.15, "lines"),
      size = text_size,
      fontface = "bold",
      seed = 42,
      max.overlaps = Inf
    )
  }

  p <- p +

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
