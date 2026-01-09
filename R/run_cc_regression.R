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

run_cc_regression <- function(obj, species = "mouse") {
  # 1. Species Check and Gene List Preparation
  species_lower <- tolower(species)

  if (species_lower == "mouse") {
    # Convert human genes to Title Case for Mouse (e.g., MKI67 -> Mki67)
    s_genes <- stringr::str_to_title(Seurat::cc.genes.updated.2019$s.genes)
    g2m_genes <- stringr::str_to_title(Seurat::cc.genes.updated.2019$g2m.genes)
    message(paste0(">>> [run_cc_regression] Processing with mouse gene symbols."))
  } else if (species_lower == "human") {
    # Use default human genes (All Caps)
    s_genes <- Seurat::cc.genes.updated.2019$s.genes
    g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes
    message(paste0(">>> [run_cc_regression] Processing with human gene symbols."))
  } else {
    # Matches your request: Stop if species is not supported
    stop(paste0("Error: Unsupported species '", species, "'."))
  }

  # 2. Verify Gene Overlap
  # Check if genes exist in the dataset to prevent errors in CellCycleScoring
  all_cc_genes <- c(s_genes, g2m_genes)
  check_overlap <- intersect(rownames(obj), all_cc_genes)

  if (length(check_overlap) == 0) {
    stop(
      "Error: No cell cycle genes found in the Seurat object!\n",
      "Please check if rownames(obj) are Gene Symbols (e.g., Mki67) or Ensembl IDs."
    )
  } else {
    message(paste0(">>> Found ", length(check_overlap), " cell cycle genes. Proceeding to scoring..."))
  }

  # 3. Cell Cycle Scoring
  # Ensure we are using the RNA assay for scoring (raw counts)
  Seurat::DefaultAssay(obj) <- "RNA"

  obj <- Seurat::CellCycleScoring(obj,
    s.features = s_genes,
    g2m.features = g2m_genes,
    set.ident = TRUE
  )

  # 4. SCTransform Regression
  # Regress out the calculated scores
  message(">>> Running SCTransform to regress out S.Score and G2M.Score...")

  obj <- Seurat::SCTransform(obj,
    vars.to.regress = c("S.Score", "G2M.Score"),
    vst.flavor = "v2",
    verbose = FALSE
  )

  message(">>> Finished: Cell cycle has been regressed out based on SCT.
")
  return(obj)
}
