The `reduceKappa_wrapper` reduces redundant pathway enrichment results
by clustering similar categories based on their shared genes using kappa
scores. It selects representative terms for each cluster and retains
relevant pathway information. The function requires information from
four columns to perform the full clustering process: 1) unique geneset
ID, 2) descriptive name of the geneset, 3) genes returned with that
category, and 4) the significance of the category. Check the default
arguments and change accordingly for your dataframe.

There are various ways to use the function to return full or minimal
information for further interpretation. To follow this example, read in
the provided file. This file contains two trials of pathway enrichment,
with the `data_label` column indicating the trial.

    pathway_res = read_tsv("pathway-res-example.txt")

**Example 1: Reducing Redundant Pathways for a Single Trial** This
example processes enrichment results from a single trial (identified by
data\_label == “G”) and removes redundant pathways while keeping the
most significant category per cluster as the representative term.

    pathway_reduce = pathway_res |> 
      filter(data_label == "G", FDR < 0.05) |> 
      reduceKappa_wrapper()

Check out potentially helpful attributes attached to the data objects

    attributes(pathway_reduce) |> names()

    [1] "class"               "row.names"           "names"              
    [4] "n_og_terms"          "n_reduced_terms"     "genes_in_cluster_df"
    [7] "cluster_info"       

    attr(pathway_reduce, "cluster_info") |> head(5)

    # A tibble: 5 × 6
    # Groups:   cluster [2]
      cluster cluster_id cluster_term            Geneset.ID Description cluster_size
        <int> <chr>      <chr>                   <chr>      <chr>              <int>
    1       1 GO:0005759 mitochondrial matrix    GO:0005759 mitochondr…            1
    2       2 GO:0098798 mitochondrial protein-… GO:0098798 mitochondr…           28
    3       2 GO:0098798 mitochondrial protein-… GO:0045333 cellular r…           28
    4       2 GO:0098798 mitochondrial protein-… GO:0098800 inner mito…           28
    5       2 GO:0098798 mitochondrial protein-… GO:0009060 aerobic re…           28

    attr(pathway_reduce, "genes_in_cluster_df") |> head(5)

    # A tibble: 5 × 3
      cluster_id cluster_term         Genes.Returned    
      <chr>      <chr>                <chr>             
    1 GO:0005759 mitochondrial matrix ENSRNOG00000002840
    2 GO:0005759 mitochondrial matrix ENSRNOG00000017032
    3 GO:0005759 mitochondrial matrix ENSRNOG00000006930
    4 GO:0005759 mitochondrial matrix ENSRNOG00000024128
    5 GO:0005759 mitochondrial matrix ENSRNOG00000006375

**Example 2: Clustering Across Multiple Trials** To create consistent
clusters across multiple pathway enrichment trials, we use all
significant pathway-associated genes as input information. Example of
how to create clusters that apply to multiple pathway enrichment trials
by using all of the genes returned by significant categories as input
information.

    pathway_reduce = pathway_res |> 
      filter(FDR < 0.05) |> 
      # set group_slice to retain one row per remaining signficant cluster for each group
      reduceKappa_wrapper(group_slice = "data_label", geneset_id_col = "Geneset.ID", gene_col = "Genes.Returned", sig_col = "P.value", descrip_col = "Description")

**Example 3: Comparing pathway results across trials (retaining
insignificant pathways in cases that any trial returns them as
significant)** To compare pathway enrichment results across trials while
retaining information about non-significant pathways, follow the
following approach. It keeps all pathways that are significant in at
least one trial, even if they are non-significant in others

    pathway_reduce = pathway_res |> 
      # retain all categories that are significant in either "data_label" trial.
      (\(x) filter(x, Geneset.ID %in% (filter(x, FDR < 0.05) |> pull(Geneset.ID))))() |> 
      # create a new column that includes only the genes that return from significant categories
      mutate(sig_pathway_genes = ifelse(FDR < 0.05, Genes.Returned, NA)) |> 
      reduceKappa_wrapper(gene_col = "sig_pathway_genes") # set group_slice = "data_label" to keep the most significant category in the cluster per group
