# reducekappa

The `reduceKappa_wrapper` function reduces redundant pathway enrichment
results by clustering similar categories based on their shared genes
using kappa scores as described in the [Metascape
paper](https://www.nature.com/articles/s41467-019-09234-6). It selects
representative terms for each cluster and retains relevant pathway
information. The function requires information from four columns to
perform the full clustering process: 1) unique geneset ID, 2)
descriptive name of the geneset, 3) genes returned with that category,
and 4) the significance of the category. Check the default arguments and
change accordingly for your dataframe.

There are various ways to use the function to return full or minimal
information for further interpretation. To follow this example, use the
package data file `pathway_res`. This file contains two trials of
pathway enrichment, with the `data_label` column indicating the trial.

    # install.packages("pak")
    # pak::pak("weitzela/reducekappa")

    library(reducekappa)
    library(dplyr)
    data(pathway_res)

**Example 1: Reducing Redundant Pathways for a Single Trial**  
This example processes enrichment results from a single trial
(identified by data\_label == “G”) and removes redundant pathways while
keeping the most significant category per cluster as the representative
term. If you want to retain all results, do not include the
filter\_representative argument in the function call.

    pathway_reduce = pathway_res |> 
      filter(data_label == "G", FDR < 0.05) |> 
      reduceKappa_wrapper(filter_representative = TRUE)
    pathway_reduce |> 
      head(5) |> 
      mutate(across(matches("Genes.|sig_"), ~ stringr::str_trunc(.x, 15))) |> 
      knitr::kable()

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 4%" />
<col style="width: 6%" />
<col style="width: 20%" />
<col style="width: 7%" />
<col style="width: 6%" />
<col style="width: 20%" />
<col style="width: 8%" />
<col style="width: 4%" />
<col style="width: 2%" />
<col style="width: 13%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">data_label</th>
<th style="text-align: right;">cluster</th>
<th style="text-align: left;">cluster_id</th>
<th style="text-align: left;">cluster_term</th>
<th style="text-align: right;">cluster_size</th>
<th style="text-align: left;">Geneset.ID</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Genes.Returned</th>
<th style="text-align: right;">P.value</th>
<th style="text-align: right;">FDR</th>
<th style="text-align: left;">Genes.Returned_inCluster</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: right;">13</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;"><a href="GO:0045333"
class="uri">GO:0045333</a></td>
<td style="text-align: left;">cellular respiration</td>
<td style="text-align: right;">23</td>
<td style="text-align: left;"><a href="GO:0045333"
class="uri">GO:0045333</a></td>
<td style="text-align: left;">cellular respiration</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;"><a href="GO:0007005"
class="uri">GO:0007005</a></td>
<td style="text-align: left;">mitochondrion organization</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0007005"
class="uri">GO:0007005</a></td>
<td style="text-align: left;">mitochondrion organization</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">5</td>
<td style="text-align: left;"><a href="GO:0009063"
class="uri">GO:0009063</a></td>
<td style="text-align: left;">cellular amino acid catabolic process</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;"><a href="GO:0009063"
class="uri">GO:0009063</a></td>
<td style="text-align: left;">cellular amino acid catabolic process</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
</tbody>
</table>

Check out potentially helpful attributes attached to the data objects

    attributes(pathway_reduce) |> names()

    ## [1] "class"               "row.names"           "names"              
    ## [4] "n_og_terms"          "n_reduced_terms"     "genes_in_cluster_df"
    ## [7] "cluster_info"

    attr(pathway_reduce, "cluster_info") |> head(5) |> knitr::kable()

<table>
<colgroup>
<col style="width: 6%" />
<col style="width: 9%" />
<col style="width: 26%" />
<col style="width: 9%" />
<col style="width: 35%" />
<col style="width: 11%" />
</colgroup>
<thead>
<tr>
<th style="text-align: right;">cluster</th>
<th style="text-align: left;">cluster_id</th>
<th style="text-align: left;">cluster_term</th>
<th style="text-align: left;">Geneset.ID</th>
<th style="text-align: left;">Description</th>
<th style="text-align: right;">cluster_size</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: right;">3</td>
</tr>
<tr>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;"><a href="GO:0098798"
class="uri">GO:0098798</a></td>
<td style="text-align: left;">mitochondrial protein-containing
complex</td>
<td style="text-align: right;">3</td>
</tr>
<tr>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;"><a href="GO:0000313"
class="uri">GO:0000313</a></td>
<td style="text-align: left;">organellar ribosome</td>
<td style="text-align: right;">3</td>
</tr>
<tr>
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: right;">13</td>
</tr>
<tr>
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: left;"><a href="GO:0044282"
class="uri">GO:0044282</a></td>
<td style="text-align: left;">small molecule catabolic process</td>
<td style="text-align: right;">13</td>
</tr>
</tbody>
</table>

    attr(pathway_reduce, "genes_in_cluster_df") |> head(5) |> knitr::kable()

<table>
<thead>
<tr>
<th style="text-align: left;">cluster_id</th>
<th style="text-align: left;">cluster_term</th>
<th style="text-align: left;">Genes.Returned</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000002840</td>
</tr>
<tr>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000017032</td>
</tr>
<tr>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000006930</td>
</tr>
<tr>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000024128</td>
</tr>
<tr>
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
    pathway_reduce |> 
      slice_head(n = 3, by = "data_label") |> 
      mutate(across(matches("Genes.|sig_"), ~ stringr::str_trunc(.x, 15))) |> 
      knitr::kable()

<table>
<colgroup>
<col style="width: 5%" />
<col style="width: 3%" />
<col style="width: 5%" />
<col style="width: 24%" />
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 24%" />
<col style="width: 7%" />
<col style="width: 3%" />
<col style="width: 1%" />
<col style="width: 11%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">data_label</th>
<th style="text-align: right;">cluster</th>
<th style="text-align: left;">cluster_id</th>
<th style="text-align: left;">cluster_term</th>
<th style="text-align: right;">cluster_size</th>
<th style="text-align: left;">Geneset.ID</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Genes.Returned</th>
<th style="text-align: right;">P.value</th>
<th style="text-align: right;">FDR</th>
<th style="text-align: left;">Genes.Returned_inCluster</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0007606"
class="uri">GO:0007606</a></td>
<td style="text-align: left;">sensory perception of chemical
stimulus</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;"><a href="GO:0007606"
class="uri">GO:0007606</a></td>
<td style="text-align: left;">sensory perception of chemical
stimulus</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0050906"
class="uri">GO:0050906</a></td>
<td style="text-align: left;">detection of stimulus involved in sensory
perception</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0050906"
class="uri">GO:0050906</a></td>
<td style="text-align: left;">detection of stimulus involved in sensory
perception</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;"><a href="GO:0015629"
class="uri">GO:0015629</a></td>
<td style="text-align: left;">actin cytoskeleton</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;"><a href="GO:0015629"
class="uri">GO:0015629</a></td>
<td style="text-align: left;">actin cytoskeleton</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">21</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">22</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: right;">12</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">23</td>
<td style="text-align: left;"><a href="GO:0045333"
class="uri">GO:0045333</a></td>
<td style="text-align: left;">cellular respiration</td>
<td style="text-align: right;">23</td>
<td style="text-align: left;"><a href="GO:0045333"
class="uri">GO:0045333</a></td>
<td style="text-align: left;">cellular respiration</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
</tbody>
</table>

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
      reduceKappa_wrapper(gene_col = "sig_pathway_genes")

    pathway_reduce |> 
      head(5) |> 
      mutate(across(matches("Genes.|sig_"), ~ stringr::str_trunc(.x, 15))) |> 
      knitr::kable()

<table style="width:100%;">
<colgroup>
<col style="width: 5%" />
<col style="width: 3%" />
<col style="width: 5%" />
<col style="width: 25%" />
<col style="width: 6%" />
<col style="width: 5%" />
<col style="width: 25%" />
<col style="width: 7%" />
<col style="width: 3%" />
<col style="width: 1%" />
<col style="width: 8%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">data_label</th>
<th style="text-align: right;">cluster</th>
<th style="text-align: left;">cluster_id</th>
<th style="text-align: left;">cluster_term</th>
<th style="text-align: right;">cluster_size</th>
<th style="text-align: left;">Geneset.ID</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Genes.Returned</th>
<th style="text-align: right;">P.value</th>
<th style="text-align: right;">FDR</th>
<th style="text-align: left;">sig_pathway_genes</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0007606"
class="uri">GO:0007606</a></td>
<td style="text-align: left;">sensory perception of chemical
stimulus</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;"><a href="GO:0007606"
class="uri">GO:0007606</a></td>
<td style="text-align: left;">sensory perception of chemical
stimulus</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0050906"
class="uri">GO:0050906</a></td>
<td style="text-align: left;">detection of stimulus involved in sensory
perception</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0050906"
class="uri">GO:0050906</a></td>
<td style="text-align: left;">detection of stimulus involved in sensory
perception</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0007606"
class="uri">GO:0007606</a></td>
<td style="text-align: left;">sensory perception of chemical
stimulus</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;">hsa04740</td>
<td style="text-align: left;">Olfactory transduction</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0007606"
class="uri">GO:0007606</a></td>
<td style="text-align: left;">sensory perception of chemical
stimulus</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;"><a href="GO:0007608"
class="uri">GO:0007608</a></td>
<td style="text-align: left;">sensory perception of smell</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;"><a href="GO:0015629"
class="uri">GO:0015629</a></td>
<td style="text-align: left;">actin cytoskeleton</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;"><a href="GO:0015629"
class="uri">GO:0015629</a></td>
<td style="text-align: left;">actin cytoskeleton</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
</tbody>
</table>

    # set group_slice = "data_label" to keep the most significant category in the cluster per group
    pathway_reduce = pathway_res |> 
      # retain all categories that are significant in either "data_label" trial.
      (\(x) filter(x, Geneset.ID %in% (filter(x, FDR < 0.05) |> pull(Geneset.ID))))() |> 
      # create a new column that includes only the genes that return from significant categories
      mutate(sig_pathway_genes = ifelse(FDR < 0.05, Genes.Returned, NA)) |> 
      reduceKappa_wrapper(gene_col = "sig_pathway_genes", group_slice = "data_label")

    pathway_reduce |> 
      slice_head(n = 3, by = "data_label") |> 
      mutate(across(matches("Genes.|sig_"), ~ stringr::str_trunc(.x, 15))) |> 
      knitr::kable()

<table>
<colgroup>
<col style="width: 4%" />
<col style="width: 3%" />
<col style="width: 4%" />
<col style="width: 22%" />
<col style="width: 5%" />
<col style="width: 4%" />
<col style="width: 22%" />
<col style="width: 6%" />
<col style="width: 3%" />
<col style="width: 1%" />
<col style="width: 7%" />
<col style="width: 11%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;">data_label</th>
<th style="text-align: right;">cluster</th>
<th style="text-align: left;">cluster_id</th>
<th style="text-align: left;">cluster_term</th>
<th style="text-align: right;">cluster_size</th>
<th style="text-align: left;">Geneset.ID</th>
<th style="text-align: left;">Description</th>
<th style="text-align: left;">Genes.Returned</th>
<th style="text-align: right;">P.value</th>
<th style="text-align: right;">FDR</th>
<th style="text-align: left;">sig_pathway_genes</th>
<th style="text-align: left;">sig_pathway_genes_inCluster</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0007606"
class="uri">GO:0007606</a></td>
<td style="text-align: left;">sensory perception of chemical
stimulus</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;"><a href="GO:0007606"
class="uri">GO:0007606</a></td>
<td style="text-align: left;">sensory perception of chemical
stimulus</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">2</td>
<td style="text-align: left;"><a href="GO:0050906"
class="uri">GO:0050906</a></td>
<td style="text-align: left;">detection of stimulus involved in sensory
perception</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;"><a href="GO:0050906"
class="uri">GO:0050906</a></td>
<td style="text-align: left;">detection of stimulus involved in sensory
perception</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">C</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;"><a href="GO:0015629"
class="uri">GO:0015629</a></td>
<td style="text-align: left;">actin cytoskeleton</td>
<td style="text-align: right;">4</td>
<td style="text-align: left;"><a href="GO:0015629"
class="uri">GO:0015629</a></td>
<td style="text-align: left;">actin cytoskeleton</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">21</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: right;">3</td>
<td style="text-align: left;"><a href="GO:0005759"
class="uri">GO:0005759</a></td>
<td style="text-align: left;">mitochondrial matrix</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">22</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: right;">12</td>
<td style="text-align: left;"><a href="GO:0016054"
class="uri">GO:0016054</a></td>
<td style="text-align: left;">organic acid catabolic process</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
<tr>
<td style="text-align: left;">G</td>
<td style="text-align: right;">23</td>
<td style="text-align: left;"><a href="GO:0045333"
class="uri">GO:0045333</a></td>
<td style="text-align: left;">cellular respiration</td>
<td style="text-align: right;">23</td>
<td style="text-align: left;"><a href="GO:0045333"
class="uri">GO:0045333</a></td>
<td style="text-align: left;">cellular respiration</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: right;">0</td>
<td style="text-align: right;">0</td>
<td style="text-align: left;">ENSRNOG00000…</td>
<td style="text-align: left;">ENSRNOG00000…</td>
</tr>
</tbody>
</table>
