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
#' @importFrom Seurat PercentageFeatureSet VlnPlot subset NormalizeData FindVariableFeatures ScaleData
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
    # --- 1. Species Check and Matching Logic ---
    species_lower <- tolower(species)

    if (species_lower == "mouse") {
        mt_pattern <- "^mt-"
    } else if (species_lower == "human") {
        mt_pattern <- "^MT-"
    } else {
        # Stop if the species is not in the predefined list
        stop(paste0("Error: Unsupported species '", species, "'. Currently supported: 'mouse', 'human'."))
    }

    # --- 2. Input Validation ---
    if (inherits(object_list, "Seurat")) {
        object_list <- list(Sample_1 = object_list)
        input_was_single <- TRUE
    } else if (!is.list(object_list)) {
        stop("Error: 'object_list' must be a Seurat object or a list of Seurat objects.")
    }

    # Ensure every element is a Seurat object before processing
    non_seurat_idx <- which(!vapply(object_list, function(x) inherits(x, "Seurat"), logical(1)))
    if (length(non_seurat_idx) > 0) {
        stop(paste0(
            "Error: All elements of 'object_list' must be Seurat objects. Invalid at positions: ",
            paste(non_seurat_idx, collapse = ", "), "."
        ))
    }

    if (is.null(names(object_list))) {
        warning("Warning: Input list has no names. Assigning default indices as names.")
        names(object_list) <- paste0("Sample_", seq_along(object_list))
    }

    # --- 3. Define Internal Processing Function ---
    .process_single <- function(sobj, name) {
        tmp_group <- "all"
        message(paste0(">>> Processing: ", name, " [Species: ", species_lower, "]"))

        # Calculate mitochondrial percentage
        # Use tryCatch to prevent errors if no mitochondrial genes are found
        try(
            {
                sobj[["percent.mt"]] <- Seurat::PercentageFeatureSet(sobj, pattern = mt_pattern)
            },
            silent = TRUE
        )

        # Abort this sample if percent.mt was not created
        if (!"percent.mt" %in% colnames(sobj[[]])) {
            warning(paste0("Skipping sample '", name, "': failed to compute percent.mt (check mitochondrial gene pattern)."))
            return(NULL)
        }

        # Create a temporary constant group to force single-group violin plots
        sobj[["all"]] <- "all"

        # Plot before filtering
        if (plot_qc) {
            p1 <- Seurat::VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = tmp_group) +
                patchwork::plot_annotation(
                    title = paste0(name, " (Pre-filter)"),
                    subtitle = paste("Species:", species)
                )
        }

        # Filter cells (Subset)
        sobj_filtered <- subset(sobj, subset = nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_mt)

        # Plot after filtering
        if (plot_qc) {
            p2 <- Seurat::VlnPlot(sobj_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = tmp_group) +
                patchwork::plot_annotation(title = paste(name, " (Post-filter)"))

            # Combine before/after into one patchwork figure
            combined_qc <- patchwork::wrap_plots(p1, p2, ncol = 1)
            print(combined_qc)
        }

        # Normalization and Variable Features
        sobj_filtered <- Seurat::NormalizeData(sobj_filtered, verbose = FALSE)
        sobj_filtered <- Seurat::FindVariableFeatures(sobj_filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

        # ScaleData
        if (scale_all_genes) {
            all_genes <- rownames(sobj_filtered)
            sobj_filtered <- Seurat::ScaleData(sobj_filtered, features = all_genes, verbose = FALSE)
        } else {
            sobj_filtered <- Seurat::ScaleData(sobj_filtered, verbose = FALSE)
        }

        # Clean up temporary grouping column
        if ("all" %in% colnames(sobj_filtered[[]])) {
            sobj_filtered[["all"]] <- NULL
        }

        message(paste0(">>> Finished: ", name, " (Cells: ", ncol(sobj_filtered), ")\n"))
        return(sobj_filtered)
    }

    # --- 4. Batch Execution ---
    processed_list <- lapply(names(object_list), function(n) {
        .process_single(object_list[[n]], n)
    })

    # Drop samples that failed percent.mt computation
    keep_idx <- !vapply(processed_list, is.null, logical(1))
    if (any(!keep_idx)) {
        warning(paste0(
            "Skipped samples due to missing percent.mt: ",
            paste(names(object_list)[!keep_idx], collapse = ", ")
        ))
    }

    processed_list <- processed_list[keep_idx]
    names(processed_list) <- names(object_list)[keep_idx]

    if (input_was_single) {
        return(processed_list[[1]])
    }

    return(processed_list)
}
