#' Detect Seurat v5 Layers and Suggest Join
#'
#' A small utility to detect whether a Seurat object's assays are in Seurat v5
#' (Assay5) format with multiple layers, and to suggest or perform layer merge
#' via `JoinLayers()`.
#'
#' @param obj A Seurat object.
#' @param assays Character vector of assays to inspect. Defaults to common assays
#'   present in the object (RNA/SCT/integrated).
#' @param join_layers Logical. If TRUE, will automatically join layers for flagged
#'   assays using `JoinLayers(obj, assays = ...)`.
#' @param verbose Logical. Print informative messages.
#'
#' @return If `join_layers = TRUE`, returns the updated Seurat object. Otherwise
#'   returns a named list describing detection results per assay.
#'
#' @examples
#' \dontrun{
#' # Basic detection
#' res <- run_v5_test(obj)
#'
#' # Auto-merge layers for RNA and SCT if multi-layer detected
#' obj <- run_v5_test(obj, assays = c("RNA", "SCT"), join_layers = TRUE)
#' }
#'
#' @export
#' @importFrom Seurat Assays DefaultAssay GetAssayData
#' @importFrom SeuratObject LayerNames JoinLayers

run_v5_test <- function(obj, assays = NULL, join_layers = FALSE, verbose = TRUE) {
    # --- 参数与 assay 选择 ---
    if (is.null(assays)) {
        try_assays <- c("RNA", "SCT", "integrated")
        assays <- intersect(Seurat::Assays(obj), try_assays)
        if (length(assays) == 0) assays <- Seurat::Assays(obj)
    }

    results <- list()

    # --- 逐 assay 检测 v5 层 ---
    for (a in assays) {
        is_assay5 <- FALSE
        layer_names <- NULL
        multi_layer <- FALSE

        # 检测 Assay5
        is_assay5 <- tryCatch(
            {
                inherits(obj[[a]], "Assay5")
            },
            error = function(e) FALSE
        )

        # 获取层名称（优先 LayerNames，兼容缺失导出时的 Layers）
        layer_names <- tryCatch(
            SeuratObject::LayerNames(obj[[a]]),
            error = function(e1) {
                tryCatch(SeuratObject::Layers(obj[[a]]), error = function(e2) NULL)
            }
        )

        # 白名单：标准层不视为切分层
        standard_layers <- c("counts", "data", "scale.data")
        split_layer_names <- setdiff(layer_names, standard_layers)
        multi_layer <- length(split_layer_names) > 0

        if (verbose) {
            if (multi_layer) {
                msg <- paste0(
                    "▶️ [run_v5_test] Assay '", a, "' has fragmented layers: ",
                    paste(split_layer_names, collapse = ", ")
                )
                message(msg)
                message(paste0(
                    "❕ [run_v5_test] Consider joining layers: JoinLayers(obj, assays = '", a, "')"
                ))
            } else if (is_assay5) {
                message(paste0(
                    "✅ [run_v5_test] Assay '", a, "' is Assay5 with standard layers (Clean)."
                ))
            } else {
                message(paste0("❕ [run_v5_test] Assay '", a, "' is not Assay5 (Legacy/v3)."))
            }
        }

        results[[a]] <- list(is_assay5 = is_assay5, multi_layer = multi_layer, layers = layer_names)

        # 合并 Layers (仅对 Assay5 且确有多层)
        if (join_layers && multi_layer && is_assay5) {
            if (verbose) message(paste0("▶️ [run_v5_test] Joining layers for assay '", a, "'..."))
            old_default <- tryCatch(Seurat::DefaultAssay(obj), error = function(e) NULL)
            on.exit(
                {
                    if (!is.null(old_default)) {
                        Seurat::DefaultAssay(obj) <- old_default
                    }
                },
                add = TRUE
            )
            Seurat::DefaultAssay(obj) <- a
            obj <- tryCatch(
                {
                    SeuratObject::JoinLayers(obj, assays = a)
                },
                error = function(e1) {
                    # 若 Seurat 未导出 JoinLayers，则跳过 fallback 并提示
                    has_export <- tryCatch(
                        {
                            "JoinLayers" %in% getNamespaceExports("Seurat")
                        },
                        error = function(e) FALSE
                    )
                    if (isTRUE(has_export)) {
                        tryCatch(
                            {
                                Seurat::JoinLayers(obj, assays = a)
                            },
                            error = function(e2) {
                                warning(paste0(
                                    "⚠️ [run_v5_test] Failed to join layers for assay '", a, "': ",
                                    conditionMessage(e2)
                                ), call. = FALSE)
                                obj
                            }
                        )
                    } else {
                        warning(paste0(
                            "⚠️ [run_v5_test] JoinLayers not exported by Seurat; skipped fallback. Original error: ",
                            conditionMessage(e1)
                        ), call. = FALSE)
                        obj
                    }
                }
            )
            # 再检层
            layer_names_after <- tryCatch(
                SeuratObject::LayerNames(obj[[a]]),
                error = function(e1) {
                    tryCatch(SeuratObject::Layers(obj[[a]]), error = function(e2) NULL)
                }
            )
            split_remaining <- setdiff(layer_names_after, standard_layers)

            if (verbose) {
                if (length(split_remaining) == 0) {
                    message(paste0("✅ [run_v5_test] Layers joined successfully for assay '", a, "'."))
                } else {
                    message(paste0(
                        "⚠️ [run_v5_test] Layers still fragmented for assay '", a, "': ",
                        paste(split_remaining, collapse = ", ")
                    ))
                }
            }
        } else if (join_layers && !is_assay5 && verbose) {
            message(paste0("❕ [run_v5_test] Assay '", a, "' is not Assay5; skipping join."))
        } else if (join_layers && !multi_layer && verbose) {
            message(paste0("❕ [run_v5_test] No fragmented layers to join for assay '", a, "."))
        }
    }

    # 返回对象或检测结果
    if (join_layers) {
        return(obj)
    } else {
        return(results)
    }
}
