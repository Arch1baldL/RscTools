#' Batch Quality Control and Preprocessing Pipeline
#'
#' This function accepts a named list of Seurat objects and performs a standard QC workflow including:
#' calculating mitochondrial percentage, filtering cells based on QC metrics, plotting QC violin plots,
#' and running standard preprocessing steps (NormalizeData, FindVariableFeatures, and ScaleData).
#'
#' @param object_list A Seurat object or a named list of Seurat objects.
#' @param species Character string indicating the species. Currently supports "mouse" (default) and "human".
#' @param max_nFeature_RNA Numeric. Filtering threshold for maximum number of genes per cell. Default is 5000.
#' @param max_percent_mt Numeric. Filtering threshold for maximum mitochondrial percentage. Default is 3.
#' @param scale_all_genes Logical. Whether to scale all genes during 'ScaleData'.
#' If TRUE (default), all genes are scaled (compute-intensive). If FALSE, only variable features are scaled.
#' @param plot_qc Logical. Whether to print QC violin plots before and after filtering. Default is TRUE.
#'
#' @importFrom Seurat PercentageFeatureSet VlnPlot NormalizeData FindVariableFeatures ScaleData
#' @importFrom patchwork plot_annotation
#' @importFrom stats setNames
#'
#' @return A list of processed Seurat objects.
#' @export

run_qc <- function(object_list,
                   species = "mouse",
                   max_nFeature_RNA = 5000,
                   max_percent_mt = 3,
                   scale_all_genes = TRUE,
                   plot_qc = TRUE) {
    input_was_single <- FALSE
    # --- 1. 物种检测与匹配规则 ---
    species_lower <- tolower(species)

    if (species_lower == "mouse") {
        mt_pattern <- "^mt-"
    } else if (species_lower == "human") {
        mt_pattern <- "^MT-"
    } else {
        # 若物种不在预设列表中则报错
        stop(paste0("Error: Unsupported species '", species, "'. Currently supported: 'mouse', 'human'."))
    }

    # --- 2. 输入校验 ---
    if (inherits(object_list, "Seurat")) {
        object_list <- list(Sample_1 = object_list)
        input_was_single <- TRUE
    } else if (!is.list(object_list)) {
        stop("Error: 'object_list' must be a Seurat object or a list of Seurat objects.")
    }

    # 确保列表中的每个元素都是 Seurat 对象
    non_seurat_idx <- which(!vapply(object_list, function(x) inherits(x, "Seurat"), logical(1)))
    if (length(non_seurat_idx) > 0) {
        stop(paste0(
            "Error: All elements of 'object_list' must be Seurat objects. Invalid at positions: ",
            paste(non_seurat_idx, collapse = ", "), "."
        ))
    }

    if (is.null(names(object_list))) {
        warning(paste0(get_icon("warning"), "[run_qc] Input list has no names. Assigning default indices as names."), call. = FALSE)
        names(object_list) <- paste0("Sample_", seq_along(object_list))
    }

    # --- 3. 定义内部处理函数 ---
    .process_single <- function(sobj, name) {
        tmp_group <- "all"
        message(paste0(get_icon("step"), "[run_qc] Processing: ", name, " [Species: ", species_lower, "]"))

        # v5 层检测：若检测到 RNA assay 为 v5 多层，提示优先合并
        if (exists("run_v5_test")) {
            v5_res <- tryCatch(run_v5_test(sobj, assays = "RNA", join_layers = FALSE, verbose = TRUE), error = function(e) NULL)
            if (!is.null(v5_res) && is.list(v5_res) && !is.null(v5_res[["RNA"]])) {
                has_v5 <- isTRUE(v5_res[["RNA"]]$is_assay5) || isTRUE(v5_res[["RNA"]]$multi_layer)
                if (has_v5) {
                    warning(paste0(get_icon("warning"), "[run_qc] Sample '", name, "' has Seurat v5 multi-layer in RNA assay. Consider: obj <- run_v5_test(obj, assays = 'RNA', join_layers = TRUE) before QC."), call. = FALSE)
                }
            }
        }

        # 计算线粒体比例
        # 使用 try 以避免在找不到线粒体基因时抛错
        try(
            {
                sobj[["percent.mt"]] <- Seurat::PercentageFeatureSet(sobj, pattern = mt_pattern)
            },
            silent = TRUE
        )

        # 若未成功计算 percent.mt，则跳过该样本
        if (!"percent.mt" %in% colnames(sobj[[]])) {
            warning(paste0(get_icon("warning"), "[run_qc] Skipping sample '", name, "': failed to compute percent.mt (check mitochondrial gene pattern)."), call. = FALSE)
            return(NULL)
        }

        # 创建临时常量分组以便绘制单组小提琴图
        sobj[["all"]] <- "all"

        # 过滤前绘图
        if (plot_qc) {
            p1 <- Seurat::VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = tmp_group) +
                patchwork::plot_annotation(
                    title = paste0(name, " (Pre-filter)"),
                    subtitle = paste("Species:", species)
                )
        }

        # 过滤细胞（Subset）
        sobj_filtered <- Seurat::subset(sobj, subset = nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_mt)

        # 过滤后绘图
        if (plot_qc) {
            p2 <- Seurat::VlnPlot(sobj_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = tmp_group) +
                patchwork::plot_annotation(title = paste(name, " (Post-filter)"))

            # 合并前后图形
            combined_qc <- patchwork::wrap_plots(p1, p2, ncol = 1)
            print(combined_qc)
        }

        # 归一化与可变基因识别
        sobj_filtered <- Seurat::NormalizeData(sobj_filtered, verbose = FALSE)
        sobj_filtered <- Seurat::FindVariableFeatures(sobj_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

        # ScaleData 缩放
        if (scale_all_genes) {
            all_genes <- rownames(sobj_filtered)
            sobj_filtered <- Seurat::ScaleData(sobj_filtered, features = all_genes, verbose = FALSE)
        } else {
            sobj_filtered <- Seurat::ScaleData(sobj_filtered, verbose = FALSE)
        }

        # 清理临时分组列
        if ("all" %in% colnames(sobj_filtered[[]])) {
            sobj_filtered[["all"]] <- NULL
        }

        message(paste0(get_icon("completed"), "[run_qc] Finished: ", name, " (Cells: ", ncol(sobj_filtered), ")"))
        return(sobj_filtered)
    }

    # --- 4. 批处理执行 ---
    processed_list <- lapply(names(object_list), function(n) {
        .process_single(object_list[[n]], n)
    })

    # 丢弃未成功计算 percent.mt 的样本
    keep_idx <- !vapply(processed_list, is.null, logical(1))
    if (any(!keep_idx)) {
        warning(paste0(
            paste0(get_icon("warning"), "[run_qc] Skipped samples due to missing percent.mt: "),
            paste(names(object_list)[!keep_idx], collapse = ", ")
        ), call. = FALSE)
    }

    processed_list <- processed_list[keep_idx]
    names(processed_list) <- names(object_list)[keep_idx]

    if (input_was_single) {
        return(processed_list[[1]])
    }

    return(processed_list)
}
