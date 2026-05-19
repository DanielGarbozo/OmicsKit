# AGENTS.md

## Scope

These instructions apply to the entire OmicsKit repository.

OmicsKit is an R package. Current documentation/data work is for branch `feat-docs`.

## Branch Rules

- Work only on branch `feat-docs`.
- Do not modify branch `dev`.
- Before making repository changes, confirm the active branch is `feat-docs`.

## Main Goal

Create real TCGA-BRCA example data, real prerendered figures, and vignettes for:

1. DEA workflow
2. Pathway analysis
3. Pathway clustering
4. Omics layers
5. Clinical modeling

## Data Integrity Rules

- Do not invent data.
- Do not create synthetic biological results.
- Use real TCGA-BRCA-derived inputs and outputs.
- Add clear `stop()` messages when expected raw files are missing.

## Package API Rules

- Do not remove exported functions.
- Do not rename exported functions.
- Do not edit files under `R/` unless explicitly requested.

## File Location Rules

- Put heavy processing scripts in `data-raw/brca/`.
- Put package-ready R objects in `data/`.
- Put path-based files such as `.bw`, `.bed`, `.tsv`, and `.rds` in `inst/extdata/brca/`.
- Put prerendered vignette figures in `vignettes/figures/` and `vignettes/figures_PA/`.

## Path Rules

- Do not use absolute local paths such as `G:/My Drive/...` in package scripts.
- Use `here::here()` or project-relative paths.

## Data-Raw Script Rules

Each script in `data-raw/brca/` must state:

- Inputs
- Outputs
- Expected file locations

Scripts must fail clearly with `stop()` when required raw files are unavailable.

## Vignette Rules

- Heavy code in vignettes must use `eval = FALSE`.
- Vignettes should display real prerendered figures.
- Do not rerun heavy processing inside vignettes by default.

## Required Comparisons

Use these comparisons where applicable:

- Tumor vs Normal
- ER_positive vs ER_negative

## Omics-Layer Example

Use this sample for omics-layer examples:

- `sample_id = "TCGA-A1-A0SH-01"`

Reason: this sample has RNA, CNV, methylation 450k, and mutation data.

Use the BRCA1 region as the main genomic region for track examples.

## Clinical Modeling Rules

- The modeling vignette must use clinical data only.
- Do not use `omics_score`.
- Do not use PAM50, RPPA clusters, methylation clusters, CN clusters, or any multiomics cluster variables in clinical modeling examples.

## Devtools Commands

Run the following only when explicitly requested:

- `devtools::document()`
- `devtools::build_vignettes()`
- `devtools::check()`
