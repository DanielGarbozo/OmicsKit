# Plot and Compare ROC Curves for Binary Classification Models

Generates publication-ready ROC curves for one or more binary
classification models evaluated on a held-out dataset. For each model
the function computes AUC with a 95% CI using the DeLong method. When
exactly two models are supplied, a DeLong paired test is performed to
compare their AUCs and the result is annotated on the plot.

Models can be supplied either as fitted `glm` objects or as named
numeric vectors of predicted probabilities, making the function flexible
for any classifier that outputs scores.

## Usage

``` r
nice_ROC(
  models,
  data = NULL,
  outcome,
  colors = c("#E63946", "#457B9D", "#2A9D8F", "#E9C46A"),
  show_ci = TRUE,
  show_delong = TRUE,
  smooth = FALSE,
  direction = "auto",
  plot_title = "ROC Curve Comparison",
  plot_subtitle = NULL,
  theme_fn = ggplot2::theme_classic,
  return_data = FALSE
)
```

## Arguments

- models:

  A *named* list. Each element is one of:

  - A fitted `glm` object (predictions generated internally via
    `predict(model, newdata = data, type = "response")`).

  - A numeric vector of predicted probabilities of length `nrow(data)`
    (values in \[0, 1\]).

  Names are used as legend labels (e.g.,
  `list("Clinical" = fit_A, "Clinical + Omics" = fit_B)`).

- data:

  A `data.frame` or `tibble`. Test / evaluation dataset. Required when
  any element of `models` is a `glm` object.

- outcome:

  Character (length 1). Name of the binary outcome column in `data`
  (must be `0`/`1` integer or numeric; `1` = event / positive class).

- colors:

  Character vector of colour codes for the ROC curves. Recycled if
  shorter than the number of models. Default: a four-colour
  colorblind-friendly palette
  (`c("#E63946", "#457B9D", "#2A9D8F", "#E9C46A")`).

- show_ci:

  Logical. Whether to display the DeLong 95% CI for each AUC in the
  legend. Default: `TRUE`.

- show_delong:

  Logical. When exactly two models are provided, annotate the plot with
  the DeLong test p-value and delta AUC. Default: `TRUE`.

- smooth:

  Logical. Apply kernel smoothing to the ROC curves (`pROC::smooth`).
  Useful for small samples; may hide real variability on large samples.
  Default: `FALSE`.

- direction:

  Character. Direction for
  [`pROC::roc`](https://rdrr.io/pkg/pROC/man/roc.html). Usually `"auto"`
  (default). Set to `"<"` or `">"` if you need to enforce a direction.

- plot_title:

  Character string. Main title of the figure. Default:
  `"ROC Curve Comparison"`.

- plot_subtitle:

  Character string or `NULL`. Subtitle. Default: `NULL`.

- theme_fn:

  A `ggplot2` theme function (without parentheses). Applied to the plot.
  Default:
  [`ggplot2::theme_classic`](https://ggplot2.tidyverse.org/reference/ggtheme.html).

- return_data:

  Logical. If `TRUE`, returns a named list containing the ggplot object
  *and* the underlying data/statistics. If `FALSE` (default), only the
  `ggplot` object is returned.

## Value

When `return_data = FALSE` (default): a `ggplot2` object.

When `return_data = TRUE`: a named list with elements:

- `plot`:

  The `ggplot2` ROC figure.

- `auc_table`:

  A `tibble` with one row per model containing: `model`, `auc`,
  `ci_lower`, `ci_upper`, `n_cases`, `n_controls`.

- `delong_test`:

  Result of
  [`pROC::roc.test`](https://rdrr.io/pkg/pROC/man/roc.test.html)
  (DeLong) if exactly two models were supplied; otherwise `NULL`.

- `roc_objects`:

  Named list of [`pROC::roc`](https://rdrr.io/pkg/pROC/man/roc.html)
  objects for further downstream analysis.

## Details

### AUC confidence intervals

CIs are computed with the DeLong method via
[`pROC::ci.auc`](https://rdrr.io/pkg/pROC/man/ci.auc.html) which
accounts for the correlation between paired measurements (same
subjects).

### DeLong test

When two models are compared on the *same* test set the observations are
paired, so the DeLong paired test is used
(`pROC::roc.test(method = "delong", paired = TRUE)`).

### Diagonal reference line

The grey dashed diagonal represents random performance (AUC = 0.50).

## See also

[`get_glm`](https://danielgarbozo.github.io/OmicsKit/reference/get_glm.md)
for fitting the models fed into `nice_ROC`.

## Examples

``` r
if (FALSE) { # \dontrun{
# -- Passing glm objects (predictions generated internally) -----------------
roc_result <- nice_ROC(
  models      = list("Clinical only"    = fit_A,
                     "Clinical + Omics" = fit_B),
  data        = test_data,
  outcome     = "stage_advanced",
  return_data = TRUE
)

roc_result$plot
roc_result$auc_table
roc_result$delong_test

# -- Passing probability vectors directly -----------------------------------
prob_A <- predict(fit_A, newdata = test_data, type = "response")
prob_B <- predict(fit_B, newdata = test_data, type = "response")

nice_ROC(
  models  = list("Clinical" = prob_A, "Clinical + Omics" = prob_B),
  data    = test_data,
  outcome = "stage_advanced"
)

# -- Single model (no comparison) -------------------------------------------
nice_ROC(
  models      = list("Clinical only" = fit_A),
  data        = test_data,
  outcome     = "stage_advanced",
  show_delong = FALSE
)
} # }
```
