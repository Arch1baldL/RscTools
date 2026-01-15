#' Find Partner Cells from Cell-Cell Interaction Nets
#'
#' This helper filters an interaction matrix for significant source->target edges, calculates a specificity score,
#' and returns the top partners per source for either weighted strength or interaction counts.
#'
#' @param object An object containing a `net` list with `weight`, `count`, and `pval` matrices.
#' @param source Character vector of source cell types to include.
#' @param target Character vector of target cell types to include.
#' @param measure Character string, either `"weight"` or `"count"`, controlling which net matrix to analyze.
#' @param top_n Integer number of top partners to return per source.
#' @param p_threshold Numeric significance cutoff applied to the `net$pval` matrix.
#'
#' @return A tibble of top partner records per source, sorted by source and value.
#' @export
#' @importFrom dplyr filter group_by mutate arrange slice_head ungroup
#' @importFrom Seurat Idents

runcc_find_partners <- function(object,
                                source = NULL,
                                target = NULL,
                                measure = c("weight", "count"),
                                top_n = 3,
                                p_threshold = 0.05) {
    message(paste0(get_icon("step"), "[runcc_find_partners] Running..."))
    measure <- match.arg(measure)


    if (measure == "weight") {
        data_mat <- as.matrix(object@net$weight)
    } else {
        data_mat <- as.matrix(object@net$count)
    }
    pval_obj <- object@net$pval
    if (length(dim(pval_obj)) == 3) {
        pval_mat <- apply(pval_obj, c(1, 2), min, na.rm = TRUE)
        if (!is.null(dimnames(pval_obj)[1])) rownames(pval_mat) <- dimnames(pval_obj)[[1]]
        if (!is.null(dimnames(pval_obj)[2])) colnames(pval_mat) <- dimnames(pval_obj)[[2]]
    } else {
        pval_mat <- as.matrix(pval_obj)
    }

    if (length(dim(data_mat)) != 2 || length(dim(pval_mat)) != 2) {
        stop(paste0(get_icon("error"), "[runcc_find_partners] net$weight/net/count/net$pval must be 2D matrices."), call. = FALSE)
    }

    # 若缺少行列名，用顺序填充，防止后续 intersect 返回空维度
    if (is.null(rownames(data_mat))) rownames(data_mat) <- seq_len(nrow(data_mat))
    if (is.null(colnames(data_mat))) colnames(data_mat) <- seq_len(ncol(data_mat))
    if (is.null(rownames(pval_mat))) rownames(pval_mat) <- seq_len(nrow(pval_mat))
    if (is.null(colnames(pval_mat))) colnames(pval_mat) <- seq_len(ncol(pval_mat))

    all_sources <- rownames(data_mat)
    all_targets <- colnames(data_mat)

    if (is.null(source) && is.null(target)) {
        if (inherits(object, "Seurat")) {
            default_idents <- as.character(Seurat::Idents(object))
            default_levels <- unique(default_idents)
        } else {
            default_levels <- rownames(object@net$weight)
        }
        source <- default_levels
        target <- default_levels
        message(paste0(get_icon("info"), "[runcc_find_partners] source/target not specified; using all idents (", paste(default_levels, collapse = ", "), ")."))
    }

    sources <- source
    targets <- target

    # 自动补集逻辑
    if (!is.null(sources) && is.null(targets)) {
        valid_source <- intersect(sources, all_sources)
        valid_target <- setdiff(all_targets, valid_source)
        message(paste0(get_icon("info"), "[runcc_find_partners] targets not specified, using complement: ", paste(valid_target, collapse = ", ")))
    } else if (is.null(sources) && !is.null(targets)) {
        valid_target <- intersect(targets, all_targets)
        valid_source <- setdiff(all_sources, valid_target)
        message(paste0(get_icon("info"), "[runcc_find_partners] sources not specified, using complement: ", paste(valid_source, collapse = ", ")))
    } else {
        valid_source <- intersect(sources, all_sources)
        valid_target <- intersect(targets, all_targets)
    }

    # 映射到 data_mat 的行/列位置
    source_idx_data <- match(valid_source, rownames(data_mat))
    target_idx_data <- match(valid_target, colnames(data_mat))
    keep_source <- !is.na(source_idx_data)
    keep_target <- !is.na(target_idx_data)
    valid_source <- valid_source[keep_source]
    valid_target <- valid_target[keep_target]
    source_idx_data <- source_idx_data[keep_source]
    target_idx_data <- target_idx_data[keep_target]

    if (length(valid_source) == 0 || length(valid_target) == 0) {
        message(paste0(get_icon("warning"), "[runcc_find_partners] No valid sources or targets found in the dataset."))
        return(NULL)
    }

    # 尝试匹配 pval 的行/列名
    source_idx_pval <- match(valid_source, rownames(pval_mat))
    target_idx_pval <- match(valid_target, colnames(pval_mat))
    if (all(is.na(source_idx_pval))) {
        ok <- source_idx_data <= nrow(pval_mat)
        if (!any(ok)) {
            stop(paste0(get_icon("error"), "[runcc_find_partners] net$pval has no matching rows for the chosen sources."), call. = FALSE)
        }
        valid_source <- valid_source[ok]
        source_idx_data <- source_idx_data[ok]
        source_idx_pval <- source_idx_data
        message(paste0(get_icon("info"), "[runcc_find_partners] net$pval rownames missing; using positional mapping on subset."))
    } else {
        keep <- !is.na(source_idx_pval)
        valid_source <- valid_source[keep]
        source_idx_data <- source_idx_data[keep]
        source_idx_pval <- source_idx_pval[keep]
    }
    if (all(is.na(target_idx_pval))) {
        ok <- target_idx_data <= ncol(pval_mat)
        if (!any(ok)) {
            stop(paste0(get_icon("error"), "[runcc_find_partners] net$pval has no matching columns for the chosen targets."), call. = FALSE)
        }
        valid_target <- valid_target[ok]
        target_idx_data <- target_idx_data[ok]
        target_idx_pval <- target_idx_data
        message(paste0(get_icon("info"), "[runcc_find_partners] net$pval colnames missing; using positional mapping on subset."))
    } else {
        keep <- !is.na(target_idx_pval)
        valid_target <- valid_target[keep]
        target_idx_data <- target_idx_data[keep]
        target_idx_pval <- target_idx_pval[keep]
    }

    if (length(valid_source) == 0 || length(valid_target) == 0) {
        message(paste0(get_icon("warning"), "[runcc_find_partners] No valid sources or targets found after aligning with net$pval."))
        return(NULL)
    }

    # 保证 valid_source/valid_target 维度为字符向量，防止只有一个元素时降为数值
    valid_source <- as.character(valid_source)
    valid_target <- as.character(valid_target)
    sub_data <- as.matrix(data_mat[source_idx_data, target_idx_data, drop = FALSE])
    sub_pval <- as.matrix(pval_mat[source_idx_pval, target_idx_pval, drop = FALSE])

    df_long <- as.data.frame(as.table(sub_data))
    colnames(df_long) <- c("Source", "Target", "Value")

    df_pval <- as.data.frame(as.table(sub_pval))
    df_long$Pval <- df_pval$Freq

    results <- df_long %>%
        dplyr::filter(.data$Pval < p_threshold) %>%
        dplyr::filter(.data$Value > 0) %>%
        dplyr::group_by(.data$Source) %>%
        dplyr::mutate(
            Metric_Type = measure,
            Mean_Source_Value = mean(.data$Value),
            Specificity_Score = .data$Value / .data$Mean_Source_Value
        ) %>%
        dplyr::arrange(dplyr::desc(.data$Value)) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(.data$Source, dplyr::desc(.data$Value))

    message(paste0(get_icon("completed"), "[runcc_find_partners] Completed."))
    return(results)
}
