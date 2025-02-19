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

**Example 1: Reducing Redundant Pathways for a Single Trial**  
This example processes enrichment results from a single trial
(identified by data\_label == “G”) and removes redundant pathways while
keeping the most significant category per cluster as the representative
term. If you want to retain all results, do not include the
filter\_representative argument in the function call.

    pathway_reduce = pathway_res |> 
      filter(data_label == "G", FDR < 0.05) |> 
      reduceKappa_wrapper(filter_representative = TRUE)

Check out potentially helpful attributes attached to the data objects

    attributes(pathway_reduce) |> names()

    [1] "class"               "row.names"           "names"              
    [4] "n_og_terms"          "n_reduced_terms"     "genes_in_cluster_df"
    [7] "cluster_info"       

    attr(pathway_reduce, "cluster_info") |> head(5) |> knitr::kable()

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 8%" />
<col style="width: 31%" />
<col style="width: 8%" />
<col style="width: 34%" />
<col style="width: 10%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: right;">cluster</th>
<th style="text-align: left;">cluster_id</th>
<th style="text-align: left;">cluster_term</th>
<th style="text-align: left;">Geneset.ID</th>
<th style="text-align: left;">Description</th>
<th style="text-align: right;">cluster_size</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: right;">1</td>
</tr>
<tr class="even">
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0098798"
class="uri">GO:0098798</a></td>
<td style="text-align: left;">mitochondrial protein-containing
complex</td>
<td style="text-align: left;"><a href="GO:0098798"
class="uri">GO:0098798</a></td>
<td style="text-align: left;">mitochondrial protein-containing
complex</td>
<td style="text-align: right;">28</td>
</tr>
<tr class="odd">
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0098798"
class="uri">GO:0098798</a></td>
<td style="text-align: left;">mitochondrial protein-containing
complex</td>
<td style="text-align: left;"><a href="GO:0045333"
class="uri">GO:0045333</a></td>
<td style="text-align: left;">cellular respiration</td>
<td style="text-align: right;">28</td>
</tr>
<tr class="even">
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0098798"
class="uri">GO:0098798</a></td>
<td style="text-align: left;">mitochondrial protein-containing
complex</td>
<td style="text-align: left;"><a href="GO:0098800"
class="uri">GO:0098800</a></td>
<td style="text-align: left;">inner mitochondrial membrane protein
complex</td>
<td style="text-align: right;">28</td>
</tr>
<tr class="odd">
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0098798"
class="uri">GO:0098798</a></td>
<td style="text-align: left;">mitochondrial protein-containing
complex</td>
<td style="text-align: left;"><a href="GO:0009060"
class="uri">GO:0009060</a></td>
<td style="text-align: left;">aerobic respiration</td>
<td style="text-align: right;">28</td>
</tr>
</tbody>
</table>

    attr(pathway_reduce, "genes_in_cluster_df") |> head(5) |> knitr::kable()

<table>
<thead>
<tr class="header">
<th style="text-align: left;">cluster_id</th>
<th style="text-align: left;">cluster_term</th>
<th style="text-align: left;">Genes.Returned</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000002840</td>
</tr>
<tr class="even">
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000017032</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000006930</td>
</tr>
<tr class="even">
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000024128</td>
</tr>
<tr class="odd">
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000006375</td>
</tr>
</tbody>
</table>

**Example 2: Clustering Across Multiple Trials**  
To create consistent clusters across multiple pathway enrichment trials,
we use all significant pathway-associated genes as input information.
Example of how to create clusters that apply to multiple pathway
enrichment trials by using all of the genes returned by significant
categories as input information.

    pathway_reduce = pathway_res |> 
      filter(FDR < 0.05) |> 
      # set group_slice to retain one row per remaining signficant cluster for each group
      reduceKappa_wrapper(group_slice = "data_label", 
                          geneset_id_col = "Geneset.ID", gene_col = "Genes.Returned", 
                          sig_col = "P.value", descrip_col = "Description")

**Example 3: Comparing pathway results across trials (retaining
insignificant pathways in cases that any trial returns them as
significant)**  
To compare pathway enrichment results across trials while retaining
information about non-significant pathways, follow the following
approach. It keeps all pathways that are significant in at least one
trial, even if they are non-significant in others

    pathway_reduce = pathway_res |> 
      # retain all categories that are significant in either "data_label" trial.
      (\(x) filter(x, Geneset.ID %in% (filter(x, FDR < 0.05) |> pull(Geneset.ID))))() |> 
      # create a new column that includes only the genes that return from significant categories
      mutate(sig_pathway_genes = ifelse(FDR < 0.05, Genes.Returned, NA)) |> 
      reduceKappa_wrapper(gene_col = "sig_pathway_genes") # set group_slice = "data_label" to keep the most significant category in the cluster per group
