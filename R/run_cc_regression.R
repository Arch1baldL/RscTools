#' Cell Cycle Scoring and SCTransform Regression
#'
#' This function performs cell cycle scoring based on species-specific gene lists
#' and uses SCTransform to regress out cell cycle scores ("S.Score" and "G2M.Score").
#'
#' @param obj A Seurat object.
#' @param species Character string. Currently supports "mouse" (default) and "human".
#' It automatically handles gene capitalization (Title Case for mouse, All Caps for human).
#'
#' @importFrom Seurat cc.genes.updated.2019 CellCycleScoring SCTransform DefaultAssay
#' @importFrom stringr str_to_title
#'
#' @return A Seurat object with a new 'SCT' assay containing regressed data.
#' @export

run_cc_regression <- function(obj, species = "mouse", force_normalize = FALSE, backup_data = TRUE) {
  # --- 1. 物种检测与基因列表准备 ---
  species_lower <- tolower(species)

  # 物种：mouse
  if (species_lower == "mouse") {
    # 将人类基因符号转为鼠类 Title Case（如 MKI67 -> Mki67）
    s_genes <- stringr::str_to_title(Seurat::cc.genes.updated.2019$s.genes)
    g2m_genes <- stringr::str_to_title(Seurat::cc.genes.updated.2019$g2m.genes)
    message(paste0("▶️ [run_cc_regression] Processing with mouse gene symbols."))

    # 物种：human
  } else if (species_lower == "human") {
    # 使用默认的人类全大写基因列表
    s_genes <- Seurat::cc.genes.updated.2019$s.genes
    g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes
    message(paste0("▶️ [run_cc_regression] Processing with human gene symbols."))

    # 不支持的物种
  } else {
    # 直接停止
    stop(paste0("❌ [run_cc_regression] Unsupported species '", species, "'."), call. = FALSE)
  }

  # --- 2. 验证基因重叠 ---
  # 检查基因是否存在，避免 CellCycleScoring 出错
  all_cc_genes <- c(s_genes, g2m_genes)
  check_overlap <- intersect(rownames(obj), all_cc_genes)

  if (length(check_overlap) == 0) {
    stop(
      "Error: No cell cycle genes found in the Seurat object!\n",
      "Please check if rownames(obj) are Gene Symbols (e.g., Mki67) or Ensembl IDs."
    )
  } else {
    message(paste0("▶️ [run_cc_regression] Found ", length(check_overlap), " cell cycle genes. Proceeding to scoring..."))
  }

  # --- 3. 细胞周期计分 ---
  # 确保使用 RNA assay（原始计数）进行计分
  Seurat::DefaultAssay(obj) <- "RNA"

  # v5 层检测：若检测到 RNA assay 为 v5 多层，则中断并提示先合并
  if (exists("run_v5_test")) {
    v5_res <- tryCatch(run_v5_test(obj, assays = "RNA", join_layers = FALSE, verbose = TRUE), error = function(e) NULL)
    if (!is.null(v5_res) && is.list(v5_res) && !is.null(v5_res[["RNA"]])) {
      has_v5 <- isTRUE(v5_res[["RNA"]]$is_assay5) || isTRUE(v5_res[["RNA"]]$multi_layer)
      if (has_v5) {
        stop(
          paste0(
            "⚠️ [run_cc_regression] Detected Seurat v5 multi-layer in RNA assay.\n",
            "Please run: \"run_v5_test()\" first!\n"
          ),
          call. = FALSE
        )
      }
    }
  }

  # 检测并在必要时对 RNA assay 运行 NormalizeData（在进行 CellCycleScoring 之前）
  counts_layer <- tryCatch(
    {
      Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
    },
    error = function(e) NULL
  )
  rna_layer <- tryCatch(
    {
      Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
    },
    error = function(e) NULL
  )
  if (is.null(counts_layer)) {
    counts_layer <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"), error = function(e) NULL)
  }
  if (is.null(rna_layer)) {
    rna_layer <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "data"), error = function(e) NULL)
  }

  needs_norm <- TRUE
  if (!is.null(rna_layer) && !is.null(counts_layer)) {
    # 如果 LayerData 返回 list（多个 layer），取第一个元素作为兼容处理
    if (is.list(rna_layer)) rna_layer <- rna_layer[[1]]
    if (is.list(counts_layer)) counts_layer <- counts_layer[[1]]

    # 对齐基因与细胞（优先用列名对齐），保证输入到 cor() 的向量长度一致
    common_genes <- intersect(rownames(rna_layer), rownames(counts_layer))
    if (length(common_genes) > 0) {
      common_cells <- intersect(colnames(rna_layer), colnames(counts_layer))
      if (length(common_cells) >= 5) {
        rna_m <- as.matrix(rna_layer[common_genes, common_cells, drop = FALSE])
        counts_m <- as.matrix(counts_layer[common_genes, common_cells, drop = FALSE])
      } else {
        # 如果无法按列名对齐，尝试按顺序对齐（谨慎），否则放弃相关性检测
        rna_m <- as.matrix(rna_layer[common_genes, , drop = FALSE])
        counts_m <- as.matrix(counts_layer[common_genes, , drop = FALSE])
        if (ncol(rna_m) != ncol(counts_m)) {
          message("⚠️ [run_cc_regression] Unable to align 'RNA' data and counts by cell names; will normalize.")
          rna_m <- counts_m <- NULL
        }
      }

      if (!is.null(rna_m) && !is.null(counts_m)) {
        # 随机采样一些基因用于相关性检测
        if (length(common_genes) <= 50) {
          samp_genes <- common_genes
        } else {
          samp_genes <- sample(common_genes, 50)
        }
        if (length(samp_genes) >= 5) {
          cors <- vapply(samp_genes, FUN.VALUE = 0.0, FUN = function(g) {
            x <- tryCatch(as.numeric(rna_m[g, ]), error = function(e) rep(NA_real_, ncol(rna_m)))
            y <- tryCatch(as.numeric(counts_m[g, ]), error = function(e) rep(NA_real_, ncol(counts_m)))
            if (length(x) != length(y) || all(is.na(x)) || all(is.na(y))) {
              return(NA_real_)
            }
            stats::cor(x, log1p(y), use = "pairwise.complete.obs", method = "spearman")
          })
          med_cor <- median(cors, na.rm = TRUE)
          if (!is.na(med_cor) && med_cor >= 0.8) {
            needs_norm <- FALSE
          }
        }
      }
    }
  }

  if (needs_norm) {
    message("▶️ [run_cc_regression] RNA 'data' layer not found or not normalized. Normalizing now...")
    if (!is.null(rna_layer)) {
      message("⚠️ [run_cc_regression] Existing RNA 'data' layer will be overwritten by NormalizeData().")
    }
    obj <- Seurat::NormalizeData(obj, assay = "RNA", verbose = FALSE)
    # 重新获取层数据以保证同步
    rna_layer <- tryCatch(Seurat::GetAssayData(obj, assay = "RNA", slot = "data"), error = function(e) NULL)
  } else {
    message("ℹ️ [run_cc_regression] RNA 'data' layer appears normalized; skipping NormalizeData().")
  }

  obj <- Seurat::CellCycleScoring(obj,
    s.features = s_genes,
    g2m.features = g2m_genes,
    assay = "RNA",
    set.ident = TRUE
  )

  # --- 4. SCTransform 回归 ---
  # 回归掉已计算的评分
  message("▶️ [run_cc_regression] Running SCTransform to regress out S.Score and G2M.Score...")

  obj <- Seurat::SCTransform(obj,
    vars.to.regress = c("S.Score", "G2M.Score"),
    vst.flavor = "v2",
    verbose = FALSE
  )

  message("✅ [run_cc_regression] Finished: Cell cycle has been regressed out based on SCT.")

  return(obj)
}
