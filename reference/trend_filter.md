# Filter genes by patient-level trend consistency.

This function filters differentially expressed genes according to the
consistency of their expression trend across paired patients. For each
gene \\g\\ and patient \\p\\, the function calculates the mean
expression in baseline samples and condition samples:
\\\bar{X}\_{N}(g,p)\\ and \\\bar{X}\_{T}(g,p)\\.

## Usage

``` r
trend_filter(
  expr,
  sampledata,
  results,
  baseline,
  conditions,
  sample.col = "sample_id",
  patient.col = "patient_id",
  group.col = "sample_type",
  gene.col = "ensembl",
  lfc.col = "log2FoldChange",
  ratio = 1.1,
  lfc.cutoff = 0,
  scale = c("linear", "log2"),
  require_complete_pairs = FALSE,
  na.rm = FALSE,
  return_removed = TRUE
)
```

## Arguments

- expr:

  Numeric matrix or data frame of expression values with genes as rows
  and sample IDs as columns. Row names must contain gene IDs.

- sampledata:

  Data frame with sample metadata.

- results:

  Data frame or named list of data frames containing differential
  expression results. Each data frame must contain a gene ID column and
  a log2 fold-change column.

- baseline:

  Character vector identifying the baseline group or groups in
  `group.col`.

- conditions:

  Character vector or named list identifying the condition group for
  each comparison in `results`. If `results` is a named list, names in
  `conditions` should match names in `results`.

- sample.col:

  Column in `sampledata` containing sample IDs. Default is
  `"sample_id"`.

- patient.col:

  Column in `sampledata` containing patient IDs. Default is
  `"patient_id"`.

- group.col:

  Column in `sampledata` containing group labels. Default is
  `"sample_type"`.

- gene.col:

  Column in `results` containing gene IDs. Default is `"ensembl"`.

- lfc.col:

  Column in `results` containing log2 fold changes. Default is
  `"log2FoldChange"`.

- ratio:

  Numeric value greater than 1. Defines the tolerated reversal
  threshold. Default is `1.1`.

- lfc.cutoff:

  Non-negative numeric value used to define UP and DOWN regulation from
  `lfc.col`. Default is `0`.

- scale:

  Expression scale. Use `"linear"` for normalized counts, CPM, TPM or
  similar linear-scale values. Use `"log2"` for log2-transformed
  expression values. Default is `"linear"`.

- require_complete_pairs:

  Logical. If `TRUE`, the function stops when unpaired patients are
  found. If `FALSE`, only patients with both baseline and condition
  samples are used. Default is `FALSE`.

- na.rm:

  Logical. Should missing expression values be removed when calculating
  patient-level means? Default is `FALSE`.

- return_removed:

  Logical. Should the output include a vector of removed genes? Default
  is `TRUE`.

## Value

A named list containing:

- One filtered data frame per comparison.

- `TrendGenes`: unique genes passing the trend consistency filter.

- `Diagnostics`: gene-level filtering diagnostics.

- `Summary`: comparison-level filtering summary.

- `RemovedGenes`: unique genes removed by the filter, if
  `return_removed = TRUE`.

## Details

For genes classified as UP-regulated at the group level, a gene is
removed if at least one paired patient shows: \\\bar{X}\_{N}(g,p) \>=
\bar{X}\_{T}(g,p) \* ratio\\.

For genes classified as DOWN-regulated at the group level, a gene is
removed if at least one paired patient shows: \\\bar{X}\_{T}(g,p) \>=
\bar{X}\_{N}(g,p) \* ratio\\.

The direction of regulation is defined from the group-level log2
fold-change. The patient-level consistency check is performed using the
expression matrix supplied in `expr`.

The multiplicative rule `ratio = 1.1` is appropriate for linear-scale
expression values. If `expr` is log2-transformed, set `scale = "log2"`;
the function will use `log2(ratio)` as an additive threshold.

## References

Requena D. et al. Nat Commun 15, 10887 (2024).

## See also

[`detect_filter()`](https://danielgarbozo.github.io/OmicsKit/reference/detect_filter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
trend_res <- trend_filter(
  expr = brca_rna_expr_tumor_normal_filtered,
  sampledata = brca_rna_metadata_tumor_normal,
  results = list(Tumor_vs_Normal = deseq_res),
  baseline = "Solid Tissue Normal",
  conditions = c(Tumor_vs_Normal = "Primary Tumor"),
  sample.col = "sampleID",
  patient.col = "patient",
  group.col = "sample_type"
)

length(trend_res$TrendGenes)
head(trend_res$Tumor_vs_Normal)
head(trend_res$Diagnostics)
trend_res$Summary
} # }
```
