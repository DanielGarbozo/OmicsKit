# List all species datasets available in a given Ensembl version.

Queries BioMart at the specified Ensembl version and returns a data
frame of all available gene datasets. Use this to:

- Find the exact `species` string to pass to
  [`get_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/get_annotations.md).

- Confirm that a species is available in a specific Ensembl version.

## Usage

``` r
list_ensembl_species(version = "current", filter = NULL)
```

## Arguments

- version:

  Ensembl release version as a string (e.g. `"112"`, `"114"`), or
  `"current"` for the latest release. Default = `"current"`.

- filter:

  Optional. A search string to filter results by dataset name or
  description (case-insensitive). Useful for quickly finding a species
  without scrolling through hundreds of rows. Default = `NULL` (no
  filter).

## Value

A data frame with columns:

- dataset:

  The dataset identifier to pass to `species` in
  [`get_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/get_annotations.md)
  (e.g. `"hsapiens_gene_ensembl"`).

- description:

  Human-readable species name and assembly (e.g.
  `"Human genes (GRCh38.p14)"`).

- version:

  Assembly version string.

## Details

Because Ensembl releases are universal (same version number for all
species), this function simply lists which species datasets exist in the
BioMart of the requested version. Species added in later releases will
not appear in older versions.

## See also

[`list_ensembl_versions()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_versions.md),
[`get_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/get_annotations.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# All species in version 112
list_ensembl_species(version = "112")

# Find zebrafish dataset in version 112
list_ensembl_species(version = "112", filter = "zebrafish")

# Find mouse dataset
list_ensembl_species(filter = "mouse")

# Find all fish species
list_ensembl_species(filter = "fish")
} # }
```
