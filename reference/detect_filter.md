# Find detectable genes across comparisons.

This function identifies genes with measurable expression levels across
samples. Detectable genes must meet two conditions: the baseMean and
their mean normalized counts in the phenotypes of interest must be
greater than a set threshold. It returns a list of detectable genes and
the comparisons in which they can be found.

## Usage

``` r
detect_filter(
  norm.counts,
  df.BvsA,
  df.CvsA = NULL,
  df.DvsA = NULL,
  cutoffs = c(50, 50, 0),
  samples.baseline,
  samples.condition1,
  samples.condition2 = NULL,
  samples.condition3 = NULL
)
```

## Arguments

- norm.counts:

  Data frame of the normalized counts with Ensembl IDs as rows and
  Sample IDs as columns.

- df.BvsA:

  Data frame comparing the first condition to the baseline.

- df.CvsA:

  Data frame comparing the second condition to the baseline (optional).

- df.DvsA:

  Data frame comparing the third condition to the baseline (optional).

- cutoffs:

  Vector containing threshold values for baseMean, mean normalized
  counts and Log2 Fold Change; respectively. Default: c(50, 50, 0).

- samples.baseline:

  Vector of Sample IDs or indexes corresponding to the baseline.

- samples.condition1:

  Vector of Sample IDs or indexes corresponding to the first condition.

- samples.condition2:

  Vector of Sample IDs or indexes corresponding to the second condition
  (optional).

- samples.condition3:

  Vector of Sample IDs or indexes corresponding to the third condition
  (optional).

## Value

A named list. Always contains:

- `$Comparison1`: Data frame of detectable genes from `df.BvsA`.

- `$DetectGenes`: Character vector of unique detectable gene IDs across
  all comparisons.

If `df.CvsA` is provided, also contains `$Comparison2`. If `df.DvsA` is
provided, also contains `$Comparison3`.

## See also

[`nice_VSB()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_VSB.md)
to plot expression of detected genes;
[brca_rna_expr_tumor_normal_filtered](https://danielgarbozo.github.io/OmicsKit/reference/brca_rna_expr_tumor_normal_filtered.md)
for an example normalized counts matrix.
[`trend_filter()`](https://danielgarbozo.github.io/OmicsKit/reference/trend_filter.md)
to find genes with significant expression trends across conditions.

## Examples

``` r
if (FALSE) { # \dontrun{
data(brca_rna_expr_tumor_normal_filtered)
data(brca_rna_dea_tumor_vs_normal)
data(brca_rna_metadata_tumor_normal)

# detect_filter requires an "ensembl" column in the results data frame
res <- brca_rna_dea_tumor_vs_normal
colnames(res)[colnames(res) == "gene_id"] <- "ensembl"
rownames(res) <- res$ensembl

# Get sample IDs per group
meta <- brca_rna_metadata_tumor_normal
samples_normal <- meta$patient[meta$sample_type == "Solid Tissue Normal"]
samples_tumor  <- meta$patient[meta$sample_type == "Primary Tumor"]

detected <- detect_filter(
  norm.counts        = as.data.frame(brca_rna_expr_tumor_normal_filtered),
  df.BvsA            = res,
  samples.baseline   = samples_normal,
  samples.condition1 = samples_tumor,
  cutoffs            = c(50, 50, 0)
)

# Number of detectable genes
length(detected$DetectGenes)

# Subset results
head(detected$Comparison1)
} # }
```
