# Correlate paired sample profiles between two omics layers.

Computes one correlation per shared sample between two normalized
matrices. Matrices are first matched by shared sample names and shared
feature names, making the function suitable for RNA-vs-protein
comparisons after mapping protein features to gene symbols.

## Usage

``` r
crossLayerCorr(
  mat_x,
  mat_y,
  method = c("spearman", "pearson"),
  top_n = 500,
  plot = TRUE
)
```

## Arguments

- mat_x:

  Numeric matrix for layer X with features as rows and samples as
  columns. Row names and column names are required.

- mat_y:

  Numeric matrix for layer Y with features as rows and samples as
  columns. Row names and column names are required.

- method:

  Correlation method. One of `"spearman"` or `"pearson"`. Default:
  `"spearman"`.

- top_n:

  Number of most variable shared features to use. If `NULL`, all shared
  features are used. Default: `500`.

- plot:

  Logical; if `TRUE`, returns a `ggplot2` bar plot in the `plot`
  element. Default: `TRUE`.

## Value

A list with:

- `correlations`: data frame with `sample` and `correlation` columns.

- `median_r`: median sample-level correlation.

- `n_shared_samples`: number of shared samples used.

- `n_shared_features`: number of shared features before top-variable
  filtering.

- `n_features_used`: number of features used for correlation.

- `plot`: `ggplot2` object or `NULL`.

## Details

By default, the function uses the top variable shared features to reduce
noise and returns both a sorted correlation table and an optional bar
plot.

## See also

[`concordanceDE()`](https://danielgarbozo.github.io/OmicsKit/reference/concordanceDE.md)
to classify genes by cross-layer differential concordance;
[`nice_ConcordanceScatter()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_ConcordanceScatter.md)
to visualize RNA-protein differential concordance.

## Examples

``` r
set.seed(1)
mat_rna <- matrix(rnorm(60), nrow = 10)
mat_protein <- mat_rna + matrix(rnorm(60, sd = 0.4), nrow = 10)
rownames(mat_rna) <- rownames(mat_protein) <- paste0("GENE", seq_len(10))
colnames(mat_rna) <- colnames(mat_protein) <- paste0("S", seq_len(6))

res <- crossLayerCorr(mat_rna, mat_protein, top_n = 8, plot = FALSE)
head(res$correlations)
#>    sample correlation sample_ordered
#> S3     S3   0.9523810             S3
#> S2     S2   0.9285714             S2
#> S6     S6   0.9285714             S6
#> S1     S1   0.9047619             S1
#> S4     S4   0.9047619             S4
#> S5     S5   0.9047619             S5
res$median_r
#> [1] 0.9166667
```
