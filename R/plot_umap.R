#' Visualize UMAP embeddings with scatter points and optional group labels,
#' using the same styling as `plot_umap_border()` but without borders.
#'
#' @param obj A Seurat object containing embeddings.
#' @param group Metadata column for coloring points.
#' @param reduction Dimensional reduction to use (default "umap").
#' @param colors Named vector for manual color mapping. Names must match
#'   `group` values (e.g., `c(Tcell = "#1f77b4", Bcell = "#ff7f0e")`).
#'   Unspecified groups will fall back to the default palette; unnamed vectors
#'   will trigger a warning.
#' @param pt_size Point size. Default 0.5.
#' @param label_size Axis label size. Default 5.
#' @param show_labels Logical; show group labels. Default TRUE.
#' @param text_size Label text size. Default 4.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData Embeddings
#' @importFrom ggplot2 ggplot aes geom_point annotate theme_void theme coord_fixed scale_color_discrete scale_color_manual margin arrow unit
#' @importFrom dplyr group_by summarise
#' @importFrom grid arrow unit
#' @importFrom rlang .data
#' @importFrom ggrepel geom_label_repel
#' @importFrom stats median

plot_umap <- function(obj,
                      group,
                      reduction = "umap",
                      colors = NULL,
                      pt_size = 0.5,
                      label_size = 5,
                      show_labels = TRUE,
                      text_size = 4) {
    message(paste0(get_icon("step"), "[plot_umap] Plotting..."))

    # show_labels 参数校验
    if (!is.logical(show_labels) && !isTRUE(show_labels) && !isFALSE(show_labels)) {
        stop(paste0(get_icon("error"), "[plot_umap] Invalid 'show_labels' argument: must be TRUE or FALSE."), call. = FALSE)
    }

    # 1. 校验与坐标提取
    if (!reduction %in% names(obj@reductions)) {
        stop(paste0(get_icon("error"), "[plot_umap] Reduction '", reduction, "' not found in Seurat object. Available reductions: ", paste(names(obj@reductions), collapse = ", ")), call. = FALSE)
    }

    coord_names <- colnames(obj[[reduction]])[1:2]

    # 元数据列校验
    meta_cols <- colnames(obj[[]])
    if (!group %in% meta_cols) {
        stop(paste0(get_icon("error"), "[plot_umap] Metadata column '", group, "' not found."), call. = FALSE)
    }

    vars <- unique(c(coord_names, group))
    df <- Seurat::FetchData(obj, vars = vars)
    umap_data <- data.frame(
        UMAP_1 = df[[coord_names[1]]],
        UMAP_2 = df[[coord_names[2]]],
        group_for_points = df[[group]],
        check.names = FALSE
    )
    if (is.factor(umap_data$group_for_points)) umap_data$group_for_points <- droplevels(umap_data$group_for_points)

    # 2. 计算坐标轴位置
    min_x <- min(umap_data$UMAP_1)
    min_y <- min(umap_data$UMAP_2)
    max_x <- max(umap_data$UMAP_1)

    arrow_len <- (max_x - min_x) * 0.15
    offset <- 1
    origin_x <- min_x - offset
    origin_y <- min_y - offset
    label_gap <- (max_x - min_x) * 0.02

    # 3. 标签位置
    if (isTRUE(show_labels)) {
        label_data <- umap_data %>%
            dplyr::group_by(.data$group_for_points) %>%
            dplyr::summarise(
                UMAP_1 = stats::median(UMAP_1),
                UMAP_2 = stats::median(UMAP_2),
                .groups = "drop"
            )
    }

    # 4. 绘图
    p <- ggplot2::ggplot(umap_data, ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)) +

        # A. 散点
        ggplot2::geom_point(ggplot2::aes(color = .data$group_for_points), size = pt_size, alpha = 0.5)

    if (isTRUE(show_labels)) {
        p <- p + ggrepel::geom_label_repel(
            data = label_data,
            ggplot2::aes(x = UMAP_1, y = UMAP_2, label = group_for_points, color = group_for_points),
            fill = "white",
            label.size = 0.6,
            label.r = ggplot2::unit(0.15, "lines"),
            size = text_size,
            fontface = "bold",
            seed = 42,
            max.overlaps = Inf,
            show.legend = FALSE
        )
    }

    if (!is.null(colors) && is.null(names(colors))) {
        warning(paste0(get_icon("warning"), "[plot_umap] 'colors' should be a named vector mapping groups to colors."), call. = FALSE)
    }

    p <- p +

        # B. 箭头坐标轴
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

        # C. 版式
        ggplot2::theme_void() +
        ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20)) +
        ggplot2::coord_fixed(clip = "off")

    if (is.null(colors)) {
        p <- p + ggplot2::scale_color_discrete()
    } else {
        p <- p + ggplot2::scale_color_manual(values = colors)
    }

    message(paste0(get_icon("completed"), "[plot_umap] Completed."))
    return(p)
}
