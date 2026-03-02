# Functions for reducing gene sets based on kappa scores, inspired by the
# Metascape pathway analysis methodology:
# https://www.nature.com/articles/s41467-019-09234-6#Sec13
# https://metascape.org/blog/?p=252

# Internal helpers ------------------------------------------------------------

.probX = function(binary_mat) {
  # input: columns = geneset, rows = 1/0 indicating presence of gene in geneset
  t1 = binary_mat |>
    (\(x) matrix(
      rep(colSums(x), ncol(x)),
      ncol = ncol(x), nrow = ncol(x),
      byrow = TRUE, dimnames = list(c(), colnames(x))
    ))()
  t2 = t(t1)
  (t1 / nrow(binary_mat)) * (t2 / nrow(binary_mat))
}

.kappaMatrix = function(mat) {
  # input: binary matrix where geneset term IDs are columns and genes are rows.
  # 1 indicates presence of gene in geneset.
  inverse_mat = 1 - mat
  p_observed_agreement = ((t(mat) %*% mat) + (t(inverse_mat) %*% inverse_mat)) / nrow(mat)
  p_chance = .probX(mat) + .probX(inverse_mat)
  (p_observed_agreement - p_chance) / (1 - p_chance)
}

.kappa2dist = function(k) {
  # Negative k values (observed overlap less than expected by chance) are scaled
  # to 0-1 before converting to a distance object (which must range 0-1).
  min_k = min(k, na.rm = TRUE)
  max_k = max(k, na.rm = TRUE)
  k = (k - min_k) / (max_k - min_k)
  as.dist(1 - k)
}

# Exported functions ----------------------------------------------------------

#' Cluster gene sets by kappa similarity
#'
#' Computes pairwise kappa similarity scores between gene sets based on shared
#' genes, then clusters them using hierarchical clustering. This implements the
#' approach described in the
#' [Metascape paper](https://www.nature.com/articles/s41467-019-09234-6).
#'
#' This is the most flexible of the two functions — it returns cluster
#' assignments as a named vector that you can use however you like in your own
#' analysis. See [reduceKappa_wrapper()] for a higher-level function that takes
#' a full enrichment results table and returns it annotated with cluster
#' information.
#'
#' @param df A two-column data frame with geneset IDs in the first column and
#'   genes or proteins in the second column, with one gene per row. Cluster
#'   similarity is calculated based on overlap in gene members across genesets.
#' @param collapse_small_clusters Logical. If `TRUE`, clusters with fewer than
#'   2 members are merged into the nearest larger cluster based on average kappa
#'   similarity. Use with caution. Default `FALSE`.
#' @param hclust_cutoff Numeric. Height at which to cut the hierarchical
#'   clustering dendrogram to define clusters. Default `0.7`, as described in
#'   the [Metascape paper](https://www.nature.com/articles/s41467-019-09234-6).
#'
#' @return A named integer vector: names are the unique values from the first
#'   column of `df` and values are integer cluster assignments from hierarchical
#'   clustering cut at `h = hclust_cutoff`.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#' data(pathway_res)
#'
#' # Build a long gene-per-row table for one trial
#' gene_long = pathway_res |>
#'   filter(data_label == "G", FDR < 0.05) |>
#'   select(Geneset.ID, Genes.Returned) |>
#'   separate_longer_delim(Genes.Returned, delim = ", ")
#'
#' clusters = reduceKappa(gene_long)
#' head(clusters)
reduceKappa = function(df, collapse_small_clusters = FALSE, hclust_cutoff = 0.7) {
  mat = df |>
    dplyr::select(1:2) |>
    `colnames<-`(c("terms_to_summarize", "supporting_info")) |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::mutate(value = 1) |>
    tidyr::pivot_wider(
      id_cols = "supporting_info",
      names_from = "terms_to_summarize",
      values_from = "value",
      values_fill = 0
    ) |>
    tibble::column_to_rownames(var = "supporting_info") |>
    as.matrix()

  k = .kappaMatrix(mat)
  hc = hclust(.kappa2dist(k), method = "average")
  clusters = cutree(hc, h = hclust_cutoff)

  if (collapse_small_clusters) {
    small_clusters = names(table(clusters)[table(clusters) < 2])
    variables_in_small_clusters = names(clusters[clusters %in% as.numeric(small_clusters)])
    new_cluster_assignments = k |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "geneset") |>
      dplyr::mutate(cluster = clusters[.data$geneset], .after = "geneset") |>
      tibble::as_tibble() |>
      dplyr::filter(!(.data$geneset %in% variables_in_small_clusters)) |>
      dplyr::select(1:2, dplyr::any_of(variables_in_small_clusters)) |>
      dplyr::group_by(.data$cluster) |>
      dplyr::summarise(dplyr::across(-"geneset", mean)) |>
      tidyr::pivot_longer(-"cluster", values_to = "avg_dist", names_to = "small_cluster_variables") |>
      dplyr::group_by(.data$small_cluster_variables) |>
      dplyr::arrange(.data$cluster) |>
      dplyr::slice_min(.data$avg_dist, n = 1, with_ties = FALSE) |>
      dplyr::pull("cluster", name = "small_cluster_variables")
    clusters[names(new_cluster_assignments)] = unname(new_cluster_assignments)
    clusters = purrr::set_names(
      dplyr::dense_rank(clusters),
      names(clusters)
    )
  }
  clusters
}


#' Reduce and annotate pathway enrichment results
#'
#' Takes an enrichment results table and returns the same table annotated with
#' cluster information. Gene sets are clustered by kappa similarity
#' (see [reduceKappa()]), and each cluster is named after the pathway with the
#' lowest p-value in that cluster. Additional summary information is attached
#' as attributes.
#'
#' @param df A data frame of pathway enrichment results. Must contain columns
#'   identified by `geneset_id_col`, `gene_col`, `sig_col`, and `descrip_col`.
#' @param group_slice Character vector of column name(s). When comparing results
#'   from multiple analyses, set this to the grouping column(s) to retain only
#'   the most significant term per cluster per group. Default `NULL`.
#' @param geneset_id_col Name of the column containing unique geneset IDs.
#'   Default `"Geneset.ID"`.
#' @param gene_col Name of the column containing genes associated with each
#'   geneset, listed as a single character string with entries separated by
#'   a space or punctuation mark (e.g. `", "` or `"/"`). Default `"Genes.Returned"`. **Note:**
#'   only genes associated with *significant* genesets should be included as
#'   supporting information to cluster pathways.
#' @param sig_col Name of the column used to select the representative term
#'   within each cluster (smallest value = most significant). Default
#'   `"P.value"`. See also `rev_sig`.
#' @param descrip_col Name of the column containing descriptive geneset names,
#'   carried through to `cluster_term`. Default `"Description"`.
#' @param delim Delimiter separating genes in `gene_col`. If `NULL` (default),
#'   the delimiter is auto-detected. Falls back to `", "` if detection fails.
#' @param rev_sig Logical. Set to `TRUE` if `sig_col` contains `-log10`
#'   transformed p-values (larger = more significant). Values are back-transformed
#'   before ranking. Default `FALSE`.
#' @param v Logical. Print progress messages. Default `FALSE`.
#' @param filter_representative Logical. If `TRUE`, retain only the single most
#'   significant row per cluster (the representative term). Default `FALSE`.
#'
#' @return A data frame based on `df` with four additional columns prepended
#'   before `geneset_id_col`:
#'   \describe{
#'     \item{`cluster`}{Integer cluster assignment.}
#'     \item{`cluster_id`}{Geneset ID of the representative (most significant)
#'       term in the cluster.}
#'     \item{`cluster_term`}{Description of the representative term.}
#'     \item{`cluster_size`}{Number of unique genesets in the cluster.}
#'   }
#'   The following attributes are attached to the returned data frame:
#'   \describe{
#'     \item{`n_og_terms`}{Number of unique genesets in the input.}
#'     \item{`n_reduced_terms`}{Number of clusters (unique representative terms).}
#'     \item{`genes_in_cluster_df`}{Data frame of all unique genes per cluster.}
#'     \item{`cluster_info`}{Data frame with detailed cluster composition.}
#'   }
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' data(pathway_res)
#'
#' # Example 1: Reduce a single trial, keep only representative terms
#' reduced = pathway_res |>
#'   filter(data_label == "G", FDR < 0.05) |>
#'   reduceKappa_wrapper(filter_representative = TRUE)
#'
#' attr(reduced, "n_og_terms")
#' attr(reduced, "n_reduced_terms")
#'
#' # Example 2: Cross-trial clustering retaining one term per group per cluster
#' reduced_multi = pathway_res |>
#'   filter(FDR < 0.05) |>
#'   reduceKappa_wrapper(group_slice = "data_label")
reduceKappa_wrapper = function(
    df,
    group_slice = NULL,
    geneset_id_col = "Geneset.ID",
    gene_col = "Genes.Returned",
    sig_col = "P.value",
    descrip_col = "Description",
    delim = NULL,
    rev_sig = FALSE,
    v = FALSE,
    filter_representative = FALSE) {

  df = dplyr::ungroup(df)
  n_og_terms = length(unique(df[[geneset_id_col]]))

  if (n_og_terms < 2) {
    message(n_og_terms, " term", ifelse(n_og_terms == 0, "s", ""), ": not enough to cluster.")
    df = df |>
      dplyr::mutate(
        cluster = NA,
        cluster_id = !!rlang::sym(geneset_id_col),
        cluster_term = !!rlang::sym(descrip_col),
        cluster_size = n_og_terms,
        .before = !!rlang::sym(geneset_id_col)
      ) |>
      `attr<-`("n_og_terms", n_og_terms) |>
      `attr<-`("n_reduced_terms", n_og_terms)
    return(df)
  }

  if (is.null(delim)) {
    delim = grep("[[:punct:]]|[[:space:]]", df[[gene_col]], value = TRUE)[1] |>
      stringr::str_remove("[:alnum:]+(?=[:punct:]|[:space:])") |>
      stringr::str_remove("(?<=[:punct:]|[:space:])[:alnum:].*$")
    if (is.na(delim)) {
      warning("could not automatically detect deliminator, setting delim = ', '")
      delim = ", "
    }
    if (v) message("detected '", delim, "' as the gene deliminator")
  }

  clusters = df |>
    dplyr::select(dplyr::all_of(c(geneset_id_col, gene_col))) |>
    tidyr::separate_longer_delim(!!rlang::sym(gene_col), delim) |>
    reduceKappa()

  if (rev_sig) {
    negLog10 = ifelse(all(df[[sig_col]] >= 0), -1, 1)
    df = df |>
      dplyr::mutate(tmp_sig = 10^(negLog10 * !!rlang::sym(sig_col)))
    sig_col = "tmp_sig"
  }

  df = df |>
    dplyr::mutate(og_idx = dplyr::row_number()) |>
    dplyr::mutate(cluster = clusters[!!rlang::sym(geneset_id_col)], .before = !!rlang::sym(geneset_id_col)) |>
    dplyr::group_by(.data$cluster) |>
    dplyr::arrange(.data$cluster, !!rlang::sym(sig_col)) |>
    dplyr::mutate(
      cluster_id = ifelse(dplyr::row_number() == 1, !!rlang::sym(geneset_id_col), NA),
      cluster_term = ifelse(dplyr::row_number() == 1, !!rlang::sym(descrip_col), NA),
      .after = "cluster"
    ) |>
    tidyr::fill("cluster_id", "cluster_term") |>
    dplyr::mutate(cluster_size = length(unique(!!rlang::sym(geneset_id_col))), .after = "cluster_term")

  genes_in_cluster = df |>
    dplyr::ungroup() |>
    dplyr::select("cluster_id", "cluster_term", dplyr::all_of(gene_col)) |>
    tidyr::separate_longer_delim(!!rlang::sym(gene_col), delim) |>
    dplyr::distinct() |>
    tidyr::drop_na()

  cluster_info = df |>
    dplyr::distinct(
      .data$cluster, .data$cluster_id, .data$cluster_term,
      !!rlang::sym(geneset_id_col), !!rlang::sym(descrip_col),
      .data$cluster_size
    )

  if (!is.null(group_slice) || filter_representative) {
    if (v) message("Filtering output dataframe to include most significant category for each group within clusters.")
    new_gene_col_nm = paste0(gene_col, "_inCluster")
    genes_in_group_cluster = df |>
      dplyr::ungroup() |>
      dplyr::select("cluster_id", "cluster_term", dplyr::all_of(group_slice), dplyr::all_of(gene_col)) |>
      tidyr::separate_longer_delim(!!rlang::sym(gene_col), delim) |>
      dplyr::distinct() |>
      tidyr::drop_na() |>
      dplyr::group_by(.data$cluster_id, .data$cluster_term, dplyr::across(dplyr::all_of(group_slice))) |>
      dplyr::summarise("{new_gene_col_nm}" := paste0(!!rlang::sym(gene_col), collapse = delim), .groups = "drop")
    df = df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c("cluster", group_slice)))) |>
      dplyr::slice_min(!!rlang::sym(sig_col), n = 1, with_ties = FALSE) |>
      dplyr::left_join(genes_in_group_cluster, by = c("cluster_id", "cluster_term", group_slice))
  }

  df = df |>
    dplyr::arrange(.data$og_idx) |>
    dplyr::select(-dplyr::any_of(c("og_idx", "tmp_sig"))) |>
    dplyr::ungroup()

  n_reduced_terms = length(unique(df[["cluster_id"]]))

  df |>
    `attr<-`("n_og_terms", n_og_terms) |>
    `attr<-`("n_reduced_terms", n_reduced_terms) |>
    `attr<-`("genes_in_cluster_df", genes_in_cluster) |>
    `attr<-`("cluster_info", cluster_info)
}
