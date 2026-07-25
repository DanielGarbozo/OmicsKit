# Plot genome-wide circos tracks for multi-omics data

Builds a circular genome plot with optional data tracks, metadata rings,
and genomic links using `circlize`.

## Usage

``` r
nice_circos(
  genome_build = "hg38",
  data_tracks = list(),
  track_types = NULL,
  metadata = NULL,
  meta_colors = NULL,
  track_colors = NULL,
  ideogram = TRUE,
  link_data = NULL,
  export_pdf = NULL,
  custom_cytoband = NULL,
  chromosome_index = NULL,
  plot_type = NULL,
  show_legend = FALSE,
  show_labels = FALSE,
  track_height = 0.1,
  metadata_height = 0.08
)
```

## Arguments

- genome_build:

  Character string. Genome build used to initialize the plot, such as
  `"hg38"`, `"hg19"`, `"mm10"`, `"ce10"`, `"ce11"`, or `"dm6"`.

- data_tracks:

  Named list of data frames, `GRanges`, file paths, or track config
  lists. Data frames must include `chr`, `start`, `end`, and optionally
  `value` and `name`.

- track_types:

  Optional character vector with one type per track. Supported values
  are `"histogram"`, `"line"`, `"points"`, `"scatter"`, and
  `"mutation"`. If `NULL`, types are inferred from track names.

- metadata:

  Optional named vector, one-row data frame, or genomic data frame with
  `chr`, `start`, `end`, and `group` columns for metadata rings.

- meta_colors:

  Optional named vector of colors for metadata groups. Names should
  correspond to values in `metadata$group`.

- track_colors:

  Optional named vector or list of colors keyed by track name.

- ideogram:

  Logical. If `TRUE`, initialize with chromosome ideogram, axis, and
  labels according to `plot_type`.

- link_data:

  Optional data frame with columns `chr1`, `start1`, `end1`, `chr2`,
  `start2`, and `end2`.

- export_pdf:

  Optional PDF file path. If `NULL`, draws on the active graphics
  device.

- custom_cytoband:

  Optional cytoband data frame with columns `chr`, `start`, `end`,
  `name`, and `gieStain`. Unnamed five-column UCSC-style cytobands are
  also accepted.

- chromosome_index:

  Optional character vector of chromosomes to keep and order.

- plot_type:

  Optional character vector passed as `plotType` to
  [`circlize::circos.initializeWithIdeogram()`](https://rdrr.io/pkg/circlize/man/circos.initializeWithIdeogram.html).
  Overrides the default from
  [`resolve_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/resolve_genome.md).

- show_legend:

  Logical. If `TRUE`, draw a simple legend using `ComplexHeatmap`.

- show_labels:

  Logical. If `TRUE`, label mutation tracks using the `name` column when
  present.

- track_height:

  Numeric scalar. Default height for data tracks.

- metadata_height:

  Numeric scalar. Height for metadata tracks.

## Value

Invisibly returns a list with:

- `tracks`:

  Normalized data tracks used for plotting.

- `metadata`:

  Normalized metadata, or `NULL`.

- `track_settings`:

  Track type, color, height, y-limits, and row count.

- `link_data`:

  Validated link data, or `NULL`.

- `genome_build`:

  Genome build used.

- `genome`:

  `circos_genome` object returned by
  [`resolve_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/resolve_genome.md),
  or `NULL` if `custom_cytoband` was used.

- `chromosomes`:

  Chromosomes shown in the plot.

- `plot_type`:

  Final plot type passed to `circlize`.

- `cytoband`:

  Cytoband used to initialize the plot.

- `pdf_opened`:

  Logical indicating whether a PDF device was opened.

## Details

Tracks must already use the selected genome build. For non-human builds
without UCSC cytobands, use
[`resolve_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/resolve_genome.md)
support for `ce10`, `ce11`, or `dm6`, or provide `custom_cytoband`.

Input tracks can be data frames, `GRanges`, file paths, or configuration
lists. File paths can point to `.bed`, `.bw`, or `.bigwig` files.
Imported `score` columns are automatically mapped to `value`.

## Track input

Each data track may be:

- a data frame with `chr`, `start`, `end`, and optionally `value`,
  `name`;

- a `GRanges` object;

- a path to `.bed`, `.bw`, or `.bigwig`;

- a config list such as
  `list(path = "x.bw", type = "line", color = "steelblue", height = 0.1)`.

## Notes

Coordinates are assumed to match `genome_build`. This function does not
perform liftOver.

For sparse regional assays, such as methylation probes or mutation calls
restricted to a single gene window, a whole-genome circos plot may make
the signal difficult to see. In that case, use `chromosome_index` to
focus on the relevant chromosome.

Mutation tracks are treated as discrete points. If no `value` column is
present, value `1` is used.

## Package requirements

This function requires `circlize`. Importing `.bed`, `.bw`, or `.bigwig`
files requires `rtracklayer`. Drawing legends requires `ComplexHeatmap`.

## References

Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize implements and
enhances circular visualization in R. Bioinformatics.
2014;30(19):2811-2812.

## See also

[`resolve_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/resolve_genome.md),
[`is_circos_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/is_circos_genome.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tracks <- list(
  RNA = list(
    path = "rna_expr_hg38.bw",
    type = "histogram",
    height = 0.12
  ),
  CNV = list(
    path = "cnv_hg38.bw",
    type = "line",
    height = 0.10
  ),
  Methylation = list(
    path = "methyl_brca1_hg38.bw",
    type = "points",
    height = 0.10
  ),
  Mutation = list(
    path = "mutations_brca1_hg38.bed",
    type = "mutation",
    height = 0.08
  )
)

nice_circos(
  genome_build = "hg38",
  data_tracks = tracks,
  chromosome_index = "chr17",
  plot_type = c("axis", "labels"),
  show_labels = TRUE,
  show_legend = TRUE,
  export_pdf = "TCGA_A1_A0SH_chr17_circos.pdf"
)
} # }

if (FALSE) { # \dontrun{
cnv <- data.frame(
  chr = c("chrI", "chrII"),
  start = c(1, 1),
  end = c(1e6, 1e6),
  value = c(0.3, -0.2)
)

nice_circos(
  genome_build = "ce10",
  data_tracks = list(CNV = cnv),
  track_types = "line"
)
} # }
```
