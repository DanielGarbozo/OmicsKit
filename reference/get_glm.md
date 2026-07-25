# Function to Fit a GLM and Return Tidy Results with Effect Sizes and FDR Correction

Fits a Generalized Linear Model (GLM) for binary, continuous, count,
proportion, rate, or positive skewed outcomes, depending on the
specified `family`. The function returns a tidy `tibble` with
coefficients, optional exponentiated effect sizes, 95% confidence
intervals, Wald test statistics, raw p-values, and
multiple-testing-adjusted q-values.

The fitted `glm` object is attached as an attribute so that it can be
passed directly to downstream functions such as
[`nice_ROC`](https://danielgarbozo.github.io/OmicsKit/reference/nice_ROC.md)
without re-fitting.

## Usage

``` r
get_glm(
  data,
  outcome,
  predictors,
  family = "binomial",
  adjust_method = "BH",
  conf_level = 0.95,
  exponentiate = NULL,
  remove_intercept = TRUE,
  verbose = FALSE
)
```

## Arguments

- data:

  A `data.frame` or `tibble` containing all outcome and predictor
  variables.

- outcome:

  Character string of length 1. Name of the outcome column. For logistic
  regression this is typically coded as `0`/`1`, although other binomial
  formats accepted by [`glm`](https://rdrr.io/r/stats/glm.html) can also
  be used.

- predictors:

  Character vector. Names of the predictor columns to include.
  Categorical variables should be `factor` or `character`; they are
  handled by [`glm`](https://rdrr.io/r/stats/glm.html) contrast coding
  automatically.

- family:

  Character string, family function, or
  [`family`](https://rdrr.io/r/stats/family.html) object passed to
  [`glm`](https://rdrr.io/r/stats/glm.html). It specifies both the
  assumed distribution of the outcome and the link function used in the
  GLM. Common options include:

  - `"binomial"` or `binomial(link = "logit")` for binary, proportion,
    or success/failure outcomes; this corresponds to logistic regression
    when the link is `"logit"`.

  - `"gaussian"` or `gaussian(link = "identity")` for approximately
    normally distributed continuous outcomes.

  - `"poisson"` or `poisson(link = "log")` for count or rate outcomes.

  - `"quasibinomial"` and `"quasipoisson"` for binomial or Poisson-type
    outcomes with overdispersion.

  - `Gamma(link = "log")` for positive, right-skewed continuous
    outcomes.

  - [`inverse.gaussian()`](https://rdrr.io/r/stats/family.html) for
    positive continuous outcomes with variance increasing strongly with
    the mean.

  Custom family objects can also be supplied. Default: `"binomial"`.

- adjust_method:

  Character string. Method for p-value adjustment passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Options: `"BH"`,
  `"bonferroni"`, `"holm"`, `"BY"`, `"fdr"`, `"none"`. Default: `"BH"`.

- conf_level:

  Numeric value in `(0, 1)`. Confidence level for the interval. Default:
  `0.95`.

- exponentiate:

  Logical or `NULL`. If `TRUE`, coefficients and confidence intervals
  are exponentiated. If `NULL`, exponentiation is applied automatically
  when the model link is `"logit"` or `"log"`. This gives odds ratios
  for logistic models with logit link, rate ratios for Poisson-type
  models with log link, and multiplicative effects for other log-link
  models. Default: `NULL`.

- remove_intercept:

  Logical. Whether to drop the `(Intercept)` row from the returned
  table. Default: `TRUE`.

- verbose:

  Logical. If `TRUE`, prints the model summary to console. Default:
  `FALSE`.

## Value

A `tibble` of class `c("get_glm_result", "tbl_df", "tbl", "data.frame")`
with the following columns:

- `term`:

  Predictor name; for factors, includes the level according to the
  contrast coding used by `glm`.

- `estimate`:

  Exponentiated coefficient if `exponentiate = TRUE`; otherwise the raw
  coefficient on the model linear predictor scale.

- `ci_lower`:

  Lower bound of the confidence interval on the same scale as
  `estimate`.

- `ci_upper`:

  Upper bound of the confidence interval on the same scale as
  `estimate`.

- `std_error`:

  Standard error of the coefficient on the linear predictor scale.

- `statistic`:

  Wald test statistic. For fixed-dispersion families such as binomial
  and Poisson this is typically a z-statistic; for families with
  estimated dispersion it may be a t-statistic.

- `p_value`:

  Two-sided p-value from the Wald test.

- `q_value`:

  P-value adjusted for multiple comparisons using the method specified
  in `adjust_method`.

- `significance`:

  Star annotation based on raw p-value: `"***"` p\<0.001, `"**"`
  p\<0.01, `"*"` p\<0.05, `"."` p\<0.10, `" "` otherwise.

The following attributes are attached to the returned object:

- `model`:

  The fitted `glm` object.

- `formula`:

  Character string of the model formula.

- `family`:

  Family name used.

- `link`:

  Link function used.

- `n_obs`:

  Number of observations used in model fitting.

- `AIC`:

  Akaike Information Criterion. May be `NA` for quasi-likelihood
  families.

- `exponentiate`:

  Logical indicating whether coefficients were exponentiated.

- `adjust_method`:

  P-value adjustment method used.

## Details

Q-values are computed across all reported terms after removing the
intercept if `remove_intercept = TRUE`. Therefore, the multiple-testing
correction pool matches the number of hypotheses shown in the returned
table.

Automatic exponentiation is based on the link function. Models with
`logit` link are reported as odds ratios; models with `log` link are
reported as multiplicative effects such as rate ratios or mean ratios.
Models with identity, probit, cloglog, inverse, or other links are not
exponentiated automatically.

## See also

[`nice_ROC`](https://danielgarbozo.github.io/OmicsKit/reference/nice_ROC.md)
for ROC curve visualisation of logistic models produced by `get_glm`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary outcome: logistic regression
res_logit <- get_glm(
  data       = train_data,
  outcome    = "stage_advanced",
  predictors = c("age", "ER", "PR", "HER2", "histology", "menopause"),
  family     = "binomial"
)

# Continuous outcome: linear model through glm
res_gaussian <- get_glm(
  data       = train_data,
  outcome    = "tumor_size",
  predictors = c("age", "ER", "PR", "HER2"),
  family     = "gaussian"
)

# Count outcome: Poisson regression
res_pois <- get_glm(
  data       = train_data,
  outcome    = "n_mutations",
  predictors = c("age", "stage", "histology"),
  family     = "poisson"
)

# Positive skewed continuous outcome
res_gamma <- get_glm(
  data       = train_data,
  outcome    = "cost",
  predictors = c("age", "stage", "treatment"),
  family     = Gamma(link = "log")
)

# Retrieve fitted model
fit <- attr(res_logit, "model")
} # }
```
