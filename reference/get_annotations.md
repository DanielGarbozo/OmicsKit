# Get gene or transcript annotations from Ensembl BioMart.

This function annotates a column of transcripts or gene IDs (ENSEMBL)
with information of the Biomart.

## Usage

``` r
get_annotations(
  ensembl_ids,
  species = "hsapiens_gene_ensembl",
  mode = "genes",
  version = "current",
  filename = "gene_annotations",
  format = "csv"
)
```

## Arguments

- ensembl_ids:

  A character vector of Ensembl gene or transcript IDs as inputs.

- species:

  The BioMart dataset identifier. Default = `"hsapiens_gene_ensembl"`
  (human). Use
  [`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md)
  to find the identifier for other organisms.

- mode:

  Either `"genes"` or `"transcripts"`. Default = `"genes"`.

- version:

  Ensembl release version as a string (e.g. `"112"`, `"114"`), or
  `"current"`. Use
  [`list_ensembl_versions()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_versions.md)
  to see available versions, and
  [`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md)
  to confirm your species exists in that version. Default = `"current"`.

- filename:

  Name of the output file without extension. If `NULL`, the result is
  returned but **not** saved to disk. Default = `"gene_annotations"`.

- format:

  Output format: `"csv"` or `"xlsx"`. Ignored when `filename = NULL`.
  Default = `"csv"`.

## Value

A data frame with one row per input ID and columns: `geneID`, `symbol`,
`biotype`, `chromosome`, `gene_start`, `gene_end`, `gene_length`,
`description`. For `mode = "transcripts"`, a `transcriptID` column is
prepended. If `filename` is not `NULL`, the data frame is also written
to disk.

## Details

The Gene information added include:

- Gene ENSEMBL ID, gene Symbol, Description, Biotype and Chromosome.

- Gene start, end and length

Annotates a vector of Ensembl gene or transcript IDs using BioMart. If
transcript IDs are provided, they are also annotated with information of
the genes to which they belong. Works for any species available in
Ensembl — use
[`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md)
to find the correct `species` value for your organism, and
[`list_ensembl_versions()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_versions.md)
to verify version availability.

**Workflow:**

1.  [`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md)
    → find dataset name.

2.  [`list_ensembl_versions()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_versions.md)
    → confirm version exists.

3.  `get_annotations(ensembl_ids, species = "drerio_gene_ensembl")`.

**Symbol resolution by species:**

- Human (`hsapiens`): uses `hgnc_symbol`, falling back to
  `external_gene_name` when HGNC symbol is absent.

- Mouse (`mmusculus`): uses `mgi_symbol`.

- Rat (`rnorvegicus`): uses `rgd_symbol`.

- Zebrafish (`drerio`): uses `zfin_id_symbol`.

- Drosophila (`dmelanogaster`): uses `flybasename_gene`.

- All other species: uses `external_gene_name` (universal fallback).

## Note

Requires an active internet connection to query the Ensembl BioMart.
`gene_length` is computed as `gene_end - gene_start + 1` (genomic span).
For TPM calculation with
[`tpm()`](https://danielgarbozo.github.io/OmicsKit/reference/tpm.md),
transcript-level lengths are more accurate.

## See also

[`list_ensembl_versions()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_versions.md),
[`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md)
to find the identifier for other organisms,
[`add_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/add_annotations.md),
[`tpm()`](https://danielgarbozo.github.io/OmicsKit/reference/tpm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Step 1 — find your species dataset
list_ensembl_species(version = "112", filter = "zebrafish")
# -> "drerio_gene_ensembl"

# Step 2 — annotate (human, default)
annotations <- get_annotations(
  ensembl_ids = c("ENSG00000141510", "ENSG00000012048"),
  species     = "hsapiens_gene_ensembl",
  version     = "112",
  filename    = NULL
)

# Mouse
annotations_mouse <- get_annotations(
  ensembl_ids = c("ENSMUSG00000059552", "ENSMUSG00000024610"),
  species     = "mmusculus_gene_ensembl",
  version     = "112",
  filename    = NULL
)

# Zebrafish
annotations_zf <- get_annotations(
  ensembl_ids = c("ENSDARG00000002333"),
  species     = "drerio_gene_ensembl",
  version     = "112",
  filename    = NULL
)
} # }
```
