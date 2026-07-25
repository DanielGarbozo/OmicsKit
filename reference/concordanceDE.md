# Classify differential results by cross-layer concordance.

Compares two differential-analysis tables from paired omics layers, such
as RNA-seq and total-protein RPPA, by merging shared genes and
classifying each gene according to statistical significance and
log-fold-change direction.

## Usage

``` r
concordanceDE(
  de_x,
  de_y,
  gene_col = "gene",
  logfc_col = "logFC",
  padj_col = "padj",
  padj_threshold = 0.05,
  logfc_threshold = 1
)
```

## Arguments

- de_x:

  Data frame for layer X. Must contain columns defined by `gene_col`,
  `logfc_col`, and `padj_col`.

- de_y:

  Data frame for layer Y. Must contain columns defined by `gene_col`,
  `logfc_col`, and `padj_col`.

- gene_col:

  Column name containing gene identifiers. Default: `"gene"`.

- logfc_col:

  Column name containing log-fold changes. Default: `"logFC"`.

- padj_col:

  Column name containing adjusted p-values/FDR values. Default:
  `"padj"`.

- padj_threshold:

  Numeric threshold for adjusted p-values. Use a single value for both
  layers or a length-two vector for layer X and layer Y, respectively.
  Default: `0.05`.

- logfc_threshold:

  Numeric absolute log-fold-change threshold. Use a single value for
  both layers or a length-two vector for layer X and layer Y,
  respectively. Default: `1`.

## Value

A list with:

- `table`: merged data frame with one `category` column.

- `concordance_score`: fraction of genes significant in both layers with
  matching log-fold-change direction.

- `summary`: table with counts per concordance category.

- `n_shared_genes`: number of genes shared by both layers.

- `n_both_significant`: number of genes significant in both layers.

## Details

Genes are assigned to one of five categories: `concordant`,
`discordant`, `only_x`, `only_y`, or `not_significant`. A concordance
score is also computed as the fraction of genes significant in both
layers that have the same effect direction.

## See also

[`nice_ConcordanceScatter()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_ConcordanceScatter.md)
to visualize the output of `concordanceDE()`;
[`crossLayerCorr()`](https://danielgarbozo.github.io/OmicsKit/reference/crossLayerCorr.md)
to estimate sample-level concordance between two omics layers.

## Examples

``` r
de_rna <- data.frame(
  gene = c("ESR1", "PGR", "EGFR", "MKI67", "GATA3"),
  logFC = c(3.2, 2.1, -1.6, -1.4, 1.8),
  padj = c(1e-6, 0.002, 0.01, 0.03, 0.04)
)

de_protein <- data.frame(
  gene = c("ESR1", "PGR", "EGFR", "MKI67", "GATA3"),
  logFC = c(0.7, 0.4, -0.5, 0.3, 0.05),
  padj = c(1e-4, 0.01, 0.02, 0.04, 0.40)
)

res <- concordanceDE(
  de_x = de_rna,
  de_y = de_protein,
  logfc_threshold = c(1, 0.2)
)
res$summary
#> 
#>      concordant      discordant          only_x          only_y not_significant 
#>               3               1               1               0               0 
res$concordance_score
#> [1] 0.75
```
