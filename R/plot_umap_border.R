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
#' @param reduction String. The dimensional reduction to use (default is "umap").
#' @param min_cells Integer. Minimum number of cells required in a group to draw a border. Groups with fewer cells will only show scatter points. Default is 10.
#' @param expand_ratio Numeric. The ratio to expand the bounding box for density calculation to avoid clipping. Default is 0.3 (30%).
#' @param pt_size Numeric. Size of the scatter points. Default is 0.5.
#' @param line_width Numeric. Width of the border lines. Default is 1.2.
#' @param label_size Numeric. Font size for UMAP axis labels. Default is 5.
#'
#' @return A ggplot object.
#' @export
#'
#' @importFrom Seurat FetchData
#' @importFrom ggplot2 ggplot aes geom_point stat_density_2d annotate theme_void theme coord_fixed scale_color_discrete margin arrow unit
#' @importFrom dplyr filter group_by group_modify ungroup bind_rows %>%
#' @importFrom grid arrow unit
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Assuming you have a Seurat object named 'seurat_obj'
#' plot_umap_border(
#'   obj = seurat_obj,
#'   group_points = "seurat_clusters",
#'   group_borders = "cell_type",
#'   min_cells = 20
#' )
#' }

plot_umap_border <- function(
    obj,
    group_points,
    group_borders,
    reduction = "umap",
    min_cells = 10,
    expand_ratio = 0.3,
    pt_size = 0.5,
    line_width = 1.2,
    label_size = 5
) {

  # 1. Validation & Coordinate Extraction
  # Check if reduction exists
  if (!reduction %in% names(obj@reductions)) {
    stop(paste0("Reduction '", reduction, "' not found in Seurat object. Available reductions: ", paste(names(obj@reductions), collapse = ", ")))
  }

  # Dynamically get the coordinate names (e.g., "umap_1", "UMAP_1", "PC_1")
  # Use [[ ]] to access the DimReduc object and get its column names
  coord_names <- colnames(obj[[reduction]])[1:2]

  # Extract UMAP coordinates and Metadata using the correct names
  umap_data <- Seurat::FetchData(obj, vars = c(coord_names, group_points, group_borders))

  # Rename columns for internal consistency (Standardize to UMAP_1/UMAP_2 for plotting logic)
  colnames(umap_data) <- c("UMAP_1", "UMAP_2", "group_for_points", "group_for_hulls")

  # 2. Filter: Keep only groups with enough cells for border drawing
  group_counts <- table(umap_data$group_for_hulls)
  valid_groups <- names(group_counts[group_counts >= min_cells])

  if (length(valid_groups) == 0) {
    stop("All groups have fewer cells than 'min_cells'. Cannot draw borders.")
  }

  umap_data_for_contour <- umap_data %>%
    dplyr::filter(.data$group_for_hulls %in% valid_groups)

  # 3. Expand boundaries (Ghost Points mechanism)
  umap_data_expanded <- umap_data_for_contour %>%
    dplyr::group_by(.data$group_for_hulls) %>%
    dplyr::group_modify(~ {
      data <- .x
      range_x <- range(data$UMAP_1)
      range_y <- range(data$UMAP_2)
      span_x <- diff(range_x) * expand_ratio
      span_y <- diff(range_y) * expand_ratio

      ghost_points <- data.frame(
        UMAP_1 = c(range_x[1]-span_x, range_x[2]+span_x, range_x[1]-span_x, range_x[2]+span_x),
        UMAP_2 = c(range_y[1]-span_y, range_y[2]+span_y, range_y[2]+span_y, range_y[1]-span_y),
        group_for_points = NA
      )
      # dplyr::bind_rows automatically handles the missing grouping column by adding it from the group context
      dplyr::bind_rows(data, ghost_points)
    }) %>%
    dplyr::ungroup()

  # 4. Calculate axis positions
  min_x <- min(umap_data$UMAP_1)
  min_y <- min(umap_data$UMAP_2)
  max_x <- max(umap_data$UMAP_1)

  arrow_len <- (max_x - min_x) * 0.15
  offset <- 1
  origin_x <- min_x - offset
  origin_y <- min_y - offset
  label_gap <- (max_x - min_x) * 0.02

  # 5. Plotting
  p <- ggplot2::ggplot(umap_data, ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)) +

    # A. Scatter points
    ggplot2::geom_point(ggplot2::aes(color = .data$group_for_points), size = pt_size, alpha = 0.5) +

    # B. Borders
    ggplot2::stat_density_2d(
      data = umap_data_expanded,
      ggplot2::aes(fill = .data$group_for_hulls, color = .data$group_for_hulls),
      geom = "polygon",
      contour_var = "ndensity",
      breaks = 0.2,
      size = line_width,
      alpha = 0.2,
      show.legend = FALSE
    ) +

    # C. Arrow Axes
    # X Axis
    ggplot2::annotate("segment", x = origin_x, xend = origin_x + arrow_len,
             y = origin_y, yend = origin_y,
             arrow = ggplot2::arrow(type = "open", length = ggplot2::unit(0.1, "inches")),
             size = 1, color = "black") +
    ggplot2::annotate("text", x = origin_x + (arrow_len * 0.05), y = origin_y - label_gap,
             label = "UMAP_1", hjust = 0, vjust = 1,
             size = label_size, fontface = "bold") +

    # Y Axis
    ggplot2::annotate("segment", x = origin_x, xend = origin_x,
             y = origin_y, yend = origin_y + arrow_len,
             arrow = ggplot2::arrow(type = "open", length = ggplot2::unit(0.1, "inches")),
             size = 1, color = "black") +
    ggplot2::annotate("text", x = origin_x - label_gap, y = origin_y + (arrow_len * 0.05),
             label = "UMAP_2", hjust = 0, vjust = 0, angle = 90,
             size = label_size, fontface = "bold") +

    # D. Aesthetics
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20)) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::scale_color_discrete(name = "Group")

  return(p)
}