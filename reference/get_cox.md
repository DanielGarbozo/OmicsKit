# Fit Cox Proportional Hazards Models and Return a Tidy Table

Fits Cox proportional hazards models against a selected survival
endpoint and returns a tidy table with hazard ratios, confidence
intervals, p-values, model metadata, and labels ready to be plotted with
[`nice_forest()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_forest.md).

## Usage

``` r
get_cox(
  data,
  time_col = "PFI.time",
  event_col = "PFI",
  vars = NULL,
  model = c("univariable", "multivariable", "adjusted"),
  adjust_vars = NULL,
  keep_adjust_terms = FALSE,
  ref_levels = NULL,
  min_n = 10,
  min_n_per_level = 5,
  min_events = 5,
  min_events_per_level = 5,
  exclude_survival_like = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  A `data.frame` containing the survival columns (`time_col`,
  `event_col`) and the variables to test.

- time_col:

  Character. Name of the column with follow-up time. The column should
  be numeric or coercible to numeric. Default: `"PFI.time"`.

- event_col:

  Character. Name of the event indicator column. It should be coded as 0
  = censored and 1 = event, or be coercible to numeric 0/1. Default:
  `"PFI"`.

- vars:

  Character vector. Variables of interest to test. If `NULL`, variables
  are selected automatically after excluding survival columns, ID-like
  columns, columns starting with `"_"`, and optionally
  survival/outcome-like columns.

- model:

  Character. Type of Cox model to fit. One of `"univariable"`,
  `"multivariable"`, or `"adjusted"`. `"univariable"` fits one model per
  variable. `"multivariable"` fits one model containing all variables in
  `vars`. `"adjusted"` fits one model per variable of interest, adjusted
  by `adjust_vars`.

- adjust_vars:

  Character vector. Covariates used when `model = "adjusted"`. Ignored
  for `"univariable"` and `"multivariable"` models. Default: `NULL`.

- keep_adjust_terms:

  Logical. If `TRUE` and `model = "adjusted"`, adjustment covariate
  terms are kept in the returned table. If `FALSE`, only terms from the
  variable of interest are returned. Default: `FALSE`.

- ref_levels:

  A named list to explicitly set reference levels for categorical
  variables. Names must be column names and values must be the desired
  reference categories. Example:
  `list(ER_Status_nature2012 = "Negative", AJCC_Stage_nature2012 = "Stage I")`.

- min_n:

  Integer. Minimum total number of complete observations required for
  each Cox model. Default: `10`.

- min_n_per_level:

  Integer. For categorical predictors, minimum number of complete
  observations required in each category. Default: `5`. Set to `0` to
  disable.

- min_events:

  Integer. Minimum total number of events required for each Cox model.
  Default: `5`. Set to `0` to disable.

- min_events_per_level:

  Integer. For categorical predictors, minimum number of events required
  in each category. This is the main safeguard against sparse categories
  and infinite Cox coefficients. Default: `5`. Set to `0` to disable.

- exclude_survival_like:

  Logical. If `TRUE` and `vars = NULL`, automatic variable selection
  excludes columns whose names suggest survival endpoints or outcome
  leakage, such as `OS`, `DSS`, `DFI`, `PFI`, `days_to`, `death`,
  `vital_status`, `followup`, or `survival`. Default: `TRUE`.

- verbose:

  Logical. If `TRUE`, prints informative messages when models or
  variables are skipped. Default: `TRUE`.

## Value

A `data.frame` with one row per model term and columns including
`model`, `model_id`, `variable`, `term`, `term_clean`, `reference`,
`HR`, `CI_low`, `CI_high`, `p.value`, `n_used`, `n_events`, and
`adjusted_for`.

## Details

The function performs several checks before fitting each model:

1.  Removes empty strings and incomplete rows.

2.  Converts character and logical predictors to factors.

3.  Drops unused factor levels.

4.  Applies explicit reference levels via `ref_levels`.

5.  Skips sparse categorical variables using `min_n_per_level` and
    `min_events_per_level`.

6.  Skips models with possible perfect separation or infinite
    coefficients.

7.  Returns exponentiated Cox coefficients as hazard ratios.

In this function, `"multivariable"` means a Cox model with multiple
predictors for one survival endpoint. This is different from
`"multivariate"`, which usually refers to multiple outcomes.

## See also

[`nice_forest`](https://danielgarbozo.github.io/OmicsKit/reference/nice_forest.md)
for plotting the tidy Cox model results returned by `get_cox()`.

[`nice_KM`](https://danielgarbozo.github.io/OmicsKit/reference/nice_KM.md)
for Kaplan-Meier survival curve visualization.

[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) and
[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html) for the underlying
Cox proportional hazards model and survival object.

[`tidy`](https://generics.r-lib.org/reference/tidy.html) for tidying
model outputs.

## Examples

``` r
if (FALSE) { # \dontrun{
# Univariable Cox models
cox_uni <- get_cox(
  data = df,
  time_col = "PFI.time",
  event_col = "PFI",
  vars = c("ER_Status_nature2012", "PAM50Call_RNAseq"),
  model = "univariable"
)

# Multivariable Cox model
cox_multi <- get_cox(
  data = df,
  time_col = "PFI.time",
  event_col = "PFI",
  vars = c("age_at_initial_pathologic_diagnosis",
           "AJCC_Stage_nature2012",
           "PAM50Call_RNAseq"),
  model = "multivariable"
)

# Adjusted Cox models
cox_adj <- get_cox(
  data = df,
  time_col = "PFI.time",
  event_col = "PFI",
  vars = c("ER_Status_nature2012", "HER2_Final_Status_nature2012"),
  adjust_vars = c("age_at_initial_pathologic_diagnosis",
                  "AJCC_Stage_nature2012"),
  model = "adjusted"
)

nice_forest(cox_adj)
} # }
```
