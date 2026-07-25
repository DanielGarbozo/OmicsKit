# Plot genomic tracks with gene annotations and optional data tracks.

Builds a multi-track genomic visualization using `Gviz`. The function
can annotate all genes in a genomic region or restrict annotations to
user-supplied gene IDs/symbols. Optional data tracks include BAM,
BigWig, and BED files.

## Usage

``` r
nice_GenomeTrack(
  region,
  genome_label = "hg38",
  organism = "hsapiens_gene_ensembl",
  ensembl_version = "current",
  annotations = NULL,
  gene_ids = NULL,
  tracks = list(),
  highlight_genes = NULL,
  gene_color = "orangered",
  other_color = "#7B68EE",
  show_transcripts = FALSE,
  track_sizes = NULL,
  export_pdf = NULL
)
```

## Arguments

- region:

  Genomic region to visualize. Either a `GRanges` object or a named
  vector with `chr`, `start`, and `end` entries (e.g.
  `c(chr = "chr1", start = 1e6, end = 2e6)`).

- genome_label:

  Assembly label used by `Gviz` (e.g. "hg38", "mm10"). This is not a
  BioMart parameter; it is only used for track metadata.

- organism:

  Ensembl BioMart dataset name (e.g. "hsapiens_gene_ensembl"). Use
  [`list_ensembl_species()`](https://danielgarbozo.github.io/OmicsKit/reference/list_ensembl_species.md)
  to find valid dataset identifiers.

- ensembl_version:

  Ensembl release version (e.g. "112") or "current". Determines which
  Ensembl archive host is queried.

- annotations:

  Optional data frame of precomputed annotations (for example from
  [`get_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/get_annotations.md)).
  Must include columns `geneID`, `symbol`, `chromosome`, `gene_start`,
  `gene_end`.

- gene_ids:

  Optional character vector of gene identifiers (Ensembl IDs or
  symbols). If provided, annotations are queried internally.

- tracks:

  Named list of file paths for data tracks. Supported formats: `bam`,
  `bw`/`bigwig`, `bed`.

- highlight_genes:

  Character vector of gene symbols to highlight.

- gene_color:

  Color for highlighted genes.

- other_color:

  Color for non-highlighted genes.

- show_transcripts:

  Logical; if `TRUE`, adds a `GeneRegionTrack` with transcripts from
  Ensembl.

- track_sizes:

  Optional numeric vector of relative heights for tracks. If provided,
  its length must match `track_list`. If `NULL`, sizes are computed
  automatically.

- export_pdf:

  Optional file path to save a PDF. If `NULL`, the plot is rendered to
  the active device.

## Value

A list of `Gviz` track objects (invisibly).

## Note

BAM files require an index file with the `.bai` suffix in the same
directory (e.g. `sample.bam.bai`). If it is missing, an error is raised.

## Examples

``` r
if (FALSE) { # \dontrun{
nice_GenomeTrack(
  region = c(chr = "chr17", start = 7e6, end = 8e6),
  genome_label = "hg38",
  organism = "hsapiens_gene_ensembl",
  ensembl_version = "current",
  tracks = list(ChIP = "chip_signal.bw")
)
} # }
```
