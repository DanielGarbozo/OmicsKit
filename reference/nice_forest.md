# Draw a Forest Plot from Cox Model Results

Creates a publication-ready forest plot from a tidy table of model
results, such as the output produced by
[`get_cox()`](https://danielgarbozo.github.io/OmicsKit/reference/get_cox.md).

## Usage

``` r
nice_forest(
  data,
  estimate_col = "HR",
  ci_low_col = "CI_low",
  ci_high_col = "CI_high",
  p_col = "p.value",
  label_col = NULL,
  variable_col = "variable",
  term_col = "term_clean",
  reference_col = "reference",
  p_display = 1,
  title = NULL,
  xlab = "Hazard Ratio (log scale)",
  log_scale = TRUE,
  vline = 1,
  color_sig = "#c0392b",
  color_ns = "#7f8c8d",
  ref_line_color = "#2c3e50",
  sort_by = c("estimate", "p.value", "input"),
  point_size = 3.5,
  base_size = 12,
  return_table = FALSE
)
```

## Arguments

- data:

  A `data.frame` containing model results. By default, it is expected to
  contain `HR`, `CI_low`, `CI_high`, and `p.value` columns.

- estimate_col:

  Character. Name of the column containing the point estimate. Default:
  `"HR"`.

- ci_low_col:

  Character. Name of the lower confidence interval column. Default:
  `"CI_low"`.

- ci_high_col:

  Character. Name of the upper confidence interval column. Default:
  `"CI_high"`.

- p_col:

  Character. Name of the p-value column. Default: `"p.value"`.

- label_col:

  Character or `NULL`. Name of a precomputed label column. If `NULL`,
  labels are built automatically from `variable_col`, `term_col`, and
  `reference_col`.

- variable_col:

  Character. Name of the variable column. Default: `"variable"`.

- term_col:

  Character. Name of the cleaned term column. Default: `"term_clean"`.

- reference_col:

  Character. Name of the reference-level column. Default: `"reference"`.

- p_display:

  Numeric. Only rows with p-value \<= `p_display` are shown. Default:
  `1` (show all rows).

- title:

  Character. Plot title. If `NULL`, a default title is used.

- xlab:

  Character. X-axis label. Default: `"Hazard Ratio (log scale)"`.

- log_scale:

  Logical. If `TRUE`, the x-axis is shown on a log10 scale. Default:
  `TRUE`.

- vline:

  Numeric. Reference line position. Default: `1`.

- color_sig:

  Color for significant points (p \< 0.05). Default: `"#c0392b"`.

- color_ns:

  Color for non-significant points. Default: `"#7f8c8d"`.

- ref_line_color:

  Color for the vertical reference line. Default: `"#2c3e50"`.

- sort_by:

  Character. How to order the plot rows. One of `"estimate"`,
  `"p.value"`, or `"input"`. Default: `"estimate"`.

- point_size:

  Numeric. Point size. Default: `3.5`.

- base_size:

  Numeric. Base font size for the ggplot theme. Default: `12`.

- return_table:

  Logical. If `TRUE`, returns a list with `$plot` and `$table`. If
  `FALSE`, returns only the plot. Default: `FALSE`.

## Value

A `ggplot` object if `return_table = FALSE`, or a named list with
`$plot` and `$table` if `return_table = TRUE`.

## Details

This function is intentionally model-agnostic. It only requires a table
with point estimates, confidence intervals, and p-values. For Cox models
generated with
[`get_cox()`](https://danielgarbozo.github.io/OmicsKit/reference/get_cox.md),
the default columns work without modification.

## See also

[`get_cox`](https://danielgarbozo.github.io/OmicsKit/reference/get_cox.md)
for fitting Cox proportional hazards models and returning tidy results
that can be passed directly to `nice_forest()`.

[`nice_KM`](https://danielgarbozo.github.io/OmicsKit/reference/nice_KM.md)
for Kaplan-Meier survival curve visualization.

[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) for the
underlying plotting system.

## Examples

``` r
if (FALSE) { # \dontrun{
cox_tab <- get_cox(
  data = df,
  time_col = "PFI.time",
  event_col = "PFI",
  vars = c("ER_Status_nature2012", "PAM50Call_RNAseq"),
  model = "univariable"
)

nice_forest(cox_tab)

nice_forest(
  cox_tab,
  p_display = 0.05,
  title = "PFI — Significant Cox Terms"
)
} # }
```
