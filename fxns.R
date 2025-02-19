# 3 functions that belong to reducing genesets based on kappa scores. I referred to https://www.nature.com/articles/s41467-019-09234-6#Sec13 and https://metascape.org/blog/?p=252 to put these functions together
.probX = function(binary_mat) {
  # input: columns = geneset, rows = 1/0 indicating presence of gene in geneset
  t1 = binary_mat |> 
    (\(x) matrix(rep(colSums(x), ncol(x)), ncol = ncol(x), nrow = ncol(x), byrow = TRUE, dimnames = list(c(), colnames(x))))()
  t2 = t(t1)
  return((t1 / nrow(binary_mat)) * (t2 / nrow(binary_mat)))
}

.kappaMatrix = function(mat) {
  # input: binary matrix where geneset term IDs are columns and genes are rows. 1 indicates presence of gene in geneset
  inverse_mat = 1 - mat
  p_observed_agreement = ((t(mat) %*% mat) + (t(inverse_mat) %*% (inverse_mat))) / nrow(mat)
  p_chance = .probX(mat) + .probX(inverse_mat)
  k = (p_observed_agreement - (p_chance)) / (1 - (p_chance))
  return(k)
}

.kappa2dist = function(k) {
  # negative k values (observed overlap is less than what is expected by chance) are addressed by scaling k to be between 0 and 1 because distance object should range from 0-1
  min_k = min(k, na.rm = TRUE)
  max_k = max(k, na.rm = TRUE)
  k = (k - min_k) / (max_k - min_k)
  k = 1 - k
  d = as.dist(k)
  return(d)
}

reduceKappa = function(df, collapse_small_clusters = FALSE) {
  # input is a dataframe where the first 2 columns are what are used for calculating similarity. the first column should be the variable you want to summarize. the second column should be the characteristics of that variable to use for summarization
  # for example, if you want to assess similarity between genesetIDs based on presence of genes, the first column should be the geneset ID and the second column should be the gene ID.
  # if you want to calculate similarity between genes based on their belonging to genesets, then the first column should be gene IDs and the second column should be geneset IDs. 
  # a named vector is returned with the geneset ID as the name and the cluster ID (number) as the value
  # collapse_small_clusters: use with caution.
  mat = df |> 
    select(1:2) |> 
    `colnames<-`(c("terms_to_summarize", "supporting_info")) |> 
    drop_na() |> 
    distinct() |> 
    mutate(value = 1) |> 
    pivot_wider(id_cols = supporting_info, names_from = "terms_to_summarize", values_from = "value", values_fill = 0) |> 
    column_to_rownames(var = "supporting_info") |> 
    as.matrix()
  k = .kappaMatrix(mat)
  hc = hclust(.kappa2dist(k), method = "average")
  clusters = cutree(hc, h = 0.7)
  
  if (collapse_small_clusters) {
    # reduce small clusters into nearest larger cluster based on closest distance to variables within the other clusters 
    small_clusters = names(table(clusters)[table(clusters)<2])
    variables_in_small_clusters = names(clusters[clusters %in% as.numeric(small_clusters)])
    new_cluster_assignments = k |>
      as.data.frame() |> rownames_to_column(var = "geneset") |>
      mutate(cluster = clusters[geneset], .after = "geneset") |>
      as_tibble() |> filter(!(geneset %in% variables_in_small_clusters)) |>
      select(1:2, any_of(variables_in_small_clusters)) |>
      group_by(cluster) |> summarise(across(-geneset, mean)) |>
      pivot_longer(-cluster, values_to = "avg_dist", names_to = "small_cluster_variables") |>
      group_by(small_cluster_variables) |>
      arrange(cluster) |>
      slice_min(avg_dist, n = 1, with_ties = FALSE) |>
      pull(cluster, name = "small_cluster_variables")
    clusters[names(new_cluster_assignments)] = unname(new_cluster_assignments)
    clusters = dense_rank(clusters) |> purrr::set_names(names(clusters))
  }
  return(clusters)
}

reduceKappa_wrapper = function(df, group_slice = NULL, geneset_id_col = "Geneset.ID", gene_col = "Genes.Returned", sig_col = "P.value", descrip_col = "Description", delim = NULL, rev_sig = FALSE, v = FALSE) {
  # Performs complete gene set reductoin, selecting representative terms and annotating clusters. Additional relevant information is provided as object attributes.
  # group_slice: IF you are comparing results from multiple analyses and you want to retain the top categories by group, then set group_slice to the column(s) that you want to keep one cluster observation per group from.
  # geneset_id_col: set this to the unique ID associated with the geneset. I prefer the GO/kegg ID over the description because they are more reliably consistent across analysis versions and there are a couple of instances of repeated description names .
  # gene_col: column of genes that will be used as information to reduce by. this assumes the genes are listed in a single character vector separated by a deliminator of some sort, which seems to be the typical output format from pathway enrichment analyses.
  # sig_col: the column that is used to determine the "representative" term of the cluster. the most significant category description/ID of the cluster is identified as the "representative" term, annotated in the cluster_id and cluster_term columns. if the input pvalues are -log10 transformed, then set rev_sig = TRUE.
  # descrip_col: column with the descriptive geneset names to be carried through to the resulting cluster_term column.
  # delim: deliminator separating genes. if not provided, function will automatically detect the deliminator. if left NULL and unable to identify a deliminator, the function will default to using ", " as the deliminator. 
  # rev_sig: see "sig_col".
  # v: print messages along the way.
  df = ungroup(df)
  n_og_terms = length(unique(df[[geneset_id_col]]))
  if (n_og_terms < 2) {
    message(n_og_terms, " term", ifelse(n_og_terms == 0, "s", ""), ": not enough to cluster.")
    df = df |> 
      mutate(cluster = NA, cluster_id = !!rlang::sym(geneset_id_col), 
             cluster_term = !!rlang::sym(descrip_col),
             cluster_size = n_og_terms, .before = !!rlang::sym(geneset_id_col)) |> 
      `attr<-`("n_og_terms", n_og_terms) |> 
      `attr<-`("n_reduced_terms", n_og_terms)
    return(df)
  }
  if (is.null(delim)) {
    delim = grep("[[:punct:]]|[[:space:]]", df[[gene_col]], value = TRUE)[1] |> 
      str_remove("[:alnum:]+(?=[:punct:]|[:space:])") |> 
      str_remove("(?<=[:punct:]|[:space:])[:alnum:].*$")
    if (is.na(delim)) {
      warning("could not automatically detect deliminator, setting delim = ', '")
      delim = ", "
    }
    if (v) message("detected '", delim, "' as the gene deliminator")
  }
  
  clusters = df |> select(all_of(c(geneset_id_col, gene_col))) |> 
    separate_longer_delim(!!rlang::sym(gene_col), delim) |> 
    reduceKappa()
  if (rev_sig) {
    # assigning names to clusters depends on the order of the significance column, from smallest to largest
    # in the case that the significance column is already log10, or -log10, then the values need to be transformed back to raw pvalues
    # if the pvalue was only log10 (all values are 0), then they do not need to be multiplied by -1
    negLog10 = ifelse(all(df[[sig_col]] >= 0), -1, 1)
    df = df |> 
      mutate(tmp_sig = 10^(negLog10 * !!rlang::sym(sig_col)))
    sig_col = "tmp_sig"
  }
  df = df |> 
    mutate(og_idx = row_number()) |> 
    mutate(cluster = clusters[!!rlang::sym(geneset_id_col)], .before = !!rlang::sym(geneset_id_col)) |> 
    group_by(cluster) |> 
    arrange(cluster, !!rlang::sym(sig_col)) |> 
    mutate(cluster_id = ifelse(row_number() == 1, !!rlang::sym(geneset_id_col), NA),
           cluster_term = ifelse(row_number() == 1, !!rlang::sym(descrip_col), NA), .after = "cluster") |> 
    fill(cluster_id, cluster_term) |>
    mutate(cluster_size = length(unique(!!rlang::sym(geneset_id_col))), .after = "cluster_term")
  genes_in_cluster = df |> 
    ungroup() |> 
    select(cluster_id, cluster_term, all_of(gene_col)) |> 
    separate_longer_delim(!!rlang::sym(gene_col), delim) |> 
    distinct() |> 
    drop_na()
  cluster_info = df |> 
    distinct(cluster, cluster_id, cluster_term, !!rlang::sym(geneset_id_col), !!rlang::sym(descrip_col), cluster_size)
  if (!is.null(group_slice)) {
    if (v) message("Filtering output dataframe to include most significant category for each group within clusters.")
    new_gene_col_nm = paste0(gene_col, "_inCluster")
    genes_in_group_cluster = df |> 
      ungroup() |> 
      select(cluster_id, cluster_term, all_of(group_slice), all_of(gene_col)) |> 
      separate_longer_delim(!!rlang::sym(gene_col), delim) |> 
      distinct() |> 
      drop_na() |> 
      group_by(cluster_id, cluster_term, across(all_of(group_slice))) |> 
      summarise("{new_gene_col_nm}" := paste0(!!rlang::sym(gene_col),  collapse = delim), .groups = "drop")
    df = df |> 
      # retain smallest pvalue per group
      group_by(across(all_of(c("cluster", group_slice)))) |> 
      slice_min(!!rlang::sym(sig_col), n = 1, with_ties = FALSE) |> 
      left_join(genes_in_group_cluster, by = c("cluster_id", "cluster_term", group_slice))
  }
  df = df |> 
    arrange(og_idx) |>
    select(-any_of(c("og_idx", "tmp_sig"))) |> 
    ungroup()
  n_reduced_terms = length(unique(df[["cluster_id"]]))
  df = df |> 
    `attr<-`("n_og_terms", n_og_terms) |> 
    `attr<-`("n_reduced_terms", n_reduced_terms) |> 
    `attr<-`("genes_in_cluster_df", genes_in_cluster) |> 
    `attr<-`("cluster_info", cluster_info)
  return(df)
}