# List available Ensembl releases and their BioMart hosts.

A convenience wrapper around
[`biomaRt::listEnsemblArchives()`](https://rdrr.io/pkg/biomaRt/man/listEnsemblArchives.html)
that returns all available Ensembl release versions with their dates and
host URLs.

## Usage

``` r
list_ensembl_versions()
```

## Value

A data frame with columns `name`, `date`, `url`, `version`, and
`current_release` (marked with `*` for the active release).

## Details

Ensembl uses a **universal release numbering system**: every release
(e.g. v112, v113) covers *all* species simultaneously. What varies
between species is whether the species is present in a given release and
which genome assembly is used. To check whether your species of interest
is available in a particular version, use
[`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md).

The recommended workflow is:

1.  Find the dataset name for your species:
    [`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md).

2.  Find available versions and confirm species presence:
    `list_ensembl_versions()` + `list_ensembl_species(version = "X")`.

3.  Annotate: `get_annotations(species = "...", version = "...")`.

## See also

[`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md),
[`get_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/get_annotations.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# See all available Ensembl versions
list_ensembl_versions()

# Then check if your species is in a specific version
list_ensembl_species(version = "112")
} # }
```
