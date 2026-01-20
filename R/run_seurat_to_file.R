#' Export Seurat counts, metadata, and reductions
#'
#' Export the counts matrix (as Matrix Market), barcodes/features, cell metadata,
#' and all available reductions from a Seurat object (v5+) into a folder for
#' downstream use (e.g., Python/Scanpy).
#'
#' @param object A Seurat object (v5+ recommended).
#' @param output_dir Output directory (created if missing).
#' @param assay Assay name to export. Default "RNA".
#' @param layer Layer name to export. Default "counts".
#' @param overwrite Logical. If FALSE (default), stops when output_dir exists and is not empty.
#'
#' @return Invisibly returns the normalized output directory path.
#' @export
#' @importFrom Seurat Assays
#' @importFrom SeuratObject JoinLayers LayerData Reductions Embeddings
#' @importFrom Matrix writeMM
#' @importFrom utils write.csv write.table

run_seurat_to_file <- function(object,
                               output_dir,
                               assay = "RNA",
                               layer = "counts",
                               overwrite = FALSE) {
    message(paste0(get_icon("step"), "[run_seurat_to_file] Running..."))

    # --- 校验 ---
    if (!inherits(object, "Seurat")) {
        stop(paste0(get_icon("error"), "[run_seurat_to_file] 'object' must be a Seurat object."), call. = FALSE)
    }
    if (missing(output_dir) || is.null(output_dir) || !nzchar(output_dir)) {
        stop(paste0(get_icon("error"), "[run_seurat_to_file] 'output_dir' cannot be empty."), call. = FALSE)
    }

    if (dir.exists(output_dir)) {
        existing <- setdiff(list.files(output_dir, all.files = TRUE, include.dirs = TRUE), c(".", ".."))
        if (length(existing) > 0 && !isTRUE(overwrite)) {
            stop(paste0(get_icon("error"), "[run_seurat_to_file] output_dir is not empty; set overwrite = TRUE to proceed."), call. = FALSE)
        }
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

    if (!assay %in% Seurat::Assays(object)) {
        stop(paste0(get_icon("error"), "[run_seurat_to_file] Assay '", assay, "' not found."), call. = FALSE)
    }

    sobj <- object
    # v5 检测与合并层（若可用 run_v5_test，则优先使用）
    if (exists("run_v5_test")) {
        v5_res <- tryCatch(run_v5_test(object, assays = assay, join_layers = FALSE, verbose = FALSE), error = function(e) NULL)
        multi_layer <- isTRUE(v5_res[[assay]]$multi_layer)
        is_assay5 <- isTRUE(v5_res[[assay]]$is_assay5)
        if (multi_layer && is_assay5) {
            message(paste0(get_icon("step"), "[run_seurat_to_file] Assay '", assay, "' has fragmented layers; joining via run_v5_test..."))
            sobj <- tryCatch(run_v5_test(object, assays = assay, join_layers = TRUE, verbose = FALSE), error = function(e) object)
        } else {
            message(paste0(get_icon("info"), "[run_seurat_to_file] Assay '", assay, "' layers clean; no join needed."))
        }
    }

    # 若仍未合并，则尝试 JoinLayers（向后兼容）
    if (identical(sobj, object)) {
        message(paste0(get_icon("step"), "[run_seurat_to_file] Joining layers for assay '", assay, "' via JoinLayers..."))
        sobj <- tryCatch(SeuratObject::JoinLayers(object, assays = assay), error = function(e) object)
    }

    counts_mat <- tryCatch(SeuratObject::LayerData(sobj, assay = assay, layer = layer), error = function(e) NULL)
    if (is.null(counts_mat) || nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
        stop(paste0(get_icon("error"), "[run_seurat_to_file] No data found in assay='", assay, "', layer='", layer, "'."), call. = FALSE)
    }

    if (!inherits(counts_mat, "dgCMatrix")) {
        counts_mat <- tryCatch(Matrix::as(counts_mat, "dgCMatrix"), error = function(e) counts_mat)
    }

    message(paste0(get_icon("step"), "[run_seurat_to_file] Writing matrix.mtx/features.tsv/barcodes.tsv..."))
    Matrix::writeMM(obj = counts_mat, file = file.path(output_dir, "matrix.mtx"))
    base::writeLines(rownames(counts_mat), con = file.path(output_dir, "features.tsv"))
    base::writeLines(colnames(counts_mat), con = file.path(output_dir, "barcodes.tsv"))

    message(paste0(get_icon("step"), "[run_seurat_to_file] Writing metadata.csv..."))
    meta_df <- sobj@meta.data
    meta_df$barcode <- rownames(meta_df)
    if (!all(meta_df$barcode == colnames(counts_mat))) {
        stop(paste0(get_icon("error"), "[run_seurat_to_file] Barcodes in metadata do not match matrix columns."), call. = FALSE)
    }
    utils::write.csv(meta_df, file = file.path(output_dir, "metadata.csv"), row.names = FALSE, quote = FALSE)

    reds <- SeuratObject::Reductions(sobj)
    if (length(reds) > 0) {
        for (red in reds) {
            message(paste0(get_icon("step"), "[run_seurat_to_file] Writing reduction: ", red, "..."))
            coords <- tryCatch(SeuratObject::Embeddings(sobj, reduction = red), error = function(e) NULL)
            if (is.null(coords) || nrow(coords) == 0 || ncol(coords) == 0) {
                warning(paste0(get_icon("warning"), "[run_seurat_to_file] Reduction '", red, "' has no embeddings; skipped."), call. = FALSE)
                next
            }
            file_name <- paste0(red, "_coords.csv")
            utils::write.csv(coords, file = file.path(output_dir, file_name), quote = FALSE)
        }
    } else {
        warning(paste0(get_icon("warning"), "[run_seurat_to_file] No reductions found; skipped."), call. = FALSE)
    }

    message(paste0(get_icon("completed"), "[run_seurat_to_file] Completed."))
    invisible(output_dir)
}
