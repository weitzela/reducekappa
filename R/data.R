#' Example pathway enrichment results
#'
#' A small subset of pathway enrichment results from two analysis trials,
#' provided for use in examples and the package vignette. The full example file
#' is available in `inst/extdata/pathway-res-example.txt`.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{`data_label`}{Character. Label identifying the analysis trial
#'     (`"G"` or `"T"`).}
#'   \item{`Geneset.ID`}{Character. Unique identifier for the gene set (e.g.
#'     GO or KEGG term ID).}
#'   \item{`Description`}{Character. Human-readable name of the gene set.}
#'   \item{`Genes.Returned`}{Character. Comma-separated list of genes
#'     associated with the gene set in this enrichment result.}
#'   \item{`P.value`}{Numeric. Nominal enrichment p-value.}
#'   \item{`FDR`}{Numeric. False discovery rate-adjusted p-value.}
#' }
#'
#' @source Subset of `pathway-res-example.txt` (available in the package
#'   repository), created via `data-raw/pathway_res.R`.
"pathway_res"
