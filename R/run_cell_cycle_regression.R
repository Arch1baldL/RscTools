#' Cell Cycle Scoring and SCTransform Regression (Force Normalize)
#'
#' Force-log-normalizes RNA before scoring, then regresses cell cycle via SCTransform.
#' Safe to run even if RNA data was already normalized.
#'
#' @param obj A Seurat object.
#' @param species Character string. Supports "mouse" (default) and "human"; mouse genes are Title Case.
#'
#' @importFrom Seurat cc.genes.updated.2019 CellCycleScoring SCTransform DefaultAssay NormalizeData
#' @importFrom stringr str_to_title
#'
#' @return A Seurat object with a new 'SCT' assay containing regressed data.
#' @export

run_cell_cycle_regression <- function(obj, species = "mouse") {
    # --- 1. 物种检测与基因列表准备 ---
    species_lower <- tolower(species)

    if (species_lower == "mouse") {
        s_genes <- stringr::str_to_title(Seurat::cc.genes.updated.2019$s.genes)
        g2m_genes <- stringr::str_to_title(Seurat::cc.genes.updated.2019$g2m.genes)
        message(sprintf("%s[run_cell_cycle_regression] Using mouse gene symbols.", get_icon("step")))
    } else if (species_lower == "human") {
        s_genes <- Seurat::cc.genes.updated.2019$s.genes
        g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes
        message(sprintf("%s[run_cell_cycle_regression] Using human gene symbols.", get_icon("step")))
    } else {
        stop(sprintf("%s[run_cell_cycle_regression] Unsupported species '%s'.", get_icon("error"), species), call. = FALSE)
    }

    # --- 2. 验证基因重叠 ---
    all_cc_genes <- c(s_genes, g2m_genes)
    check_overlap <- intersect(rownames(obj), all_cc_genes)
    if (length(check_overlap) == 0) {
        stop(
            paste0(
                sprintf("%s[run_cell_cycle_regression] No cell cycle genes found in the Seurat object.\n", get_icon("error")),
                "Please check if rownames(obj) are gene symbols (e.g., Mki67) or Ensembl IDs."
            ),
            call. = FALSE
        )
    }
    message(sprintf("%s[run_cell_cycle_regression] Found %d cell cycle genes. Proceeding...", get_icon("step"), length(check_overlap)))

    # --- 3. Assay 与层准备 ---
    Seurat::DefaultAssay(obj) <- "RNA"

    # 检测 v5 multilayers
    if (exists("run_v5_test")) {
        v5_res <- tryCatch(run_v5_test(obj, assays = "RNA", join_layers = FALSE, verbose = FALSE), error = function(e) NULL)
        if (!is.null(v5_res) && is.list(v5_res) && !is.null(v5_res[["RNA"]])) {
            has_v5_layers <- isTRUE(v5_res[["RNA"]]$multi_layer)
            if (has_v5_layers) {
                stop(
                    paste0(
                        sprintf("%s[run_cell_cycle_regression] Detected Seurat v5 fragmented layers in RNA assay.\n", get_icon("error")),
                        "Please run: run_v5_test(obj, assays = 'RNA', join_layers = TRUE) before run_cell_cycle_regression()."
                    ),
                    call. = FALSE
                )
            }
        }
    }

    # --- 4. 强制归一化 ---
    message(sprintf("%s[run_cell_cycle_regression] Normalizing RNA (LogNormalize) for scoring...", get_icon("step")))
    obj <- Seurat::NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)

    # --- 5. 细胞周期评分 ---
    message(sprintf("%s[run_cell_cycle_regression] Scoring cell cycle...", get_icon("step")))
    obj <- Seurat::CellCycleScoring(
        obj,
        s.features = s_genes,
        g2m.features = g2m_genes,
        assay = "RNA",
        set.ident = TRUE
    )

    # --- 6. SCTransform 回归 ---
    message(sprintf("%s[run_cell_cycle_regression] Running SCTransform (regressing S.Score & G2M.Score)...", get_icon("step")))
    obj <- Seurat::SCTransform(
        obj,
        vars.to.regress = c("S.Score", "G2M.Score"),
        vst.flavor = "v2",
        verbose = FALSE
    )

    message(sprintf("%s[run_cell_cycle_regression] Done.", get_icon("completed")))
    return(obj)
}
