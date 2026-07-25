# Resolve cytoband and chromosome order for circos plots

Builds a reusable genome descriptor for
[`nice_circos()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_circos.md).
For `ce10`, `ce11`, and `dm6`, chromosome sizes are provided internally,
so no UCSC download is needed. Other genome builds are resolved with
`circlize::read.cytoband(species = genome_build)`.

## Usage

``` r
resolve_genome(genome_build, chromosome_index = NULL)
```

## Arguments

- genome_build:

  Character string. Genome build name, such as `"ce10"`, `"ce11"`,
  `"dm6"`, `"hg38"`, `"hg19"`, or `"mm10"`.

- chromosome_index:

  Optional character vector of chromosomes to keep and order in the
  circos plot. If `NULL`, all available chromosomes are used.

## Value

An object of class `"circos_genome"` with elements:

- `genome_build`:

  Requested genome build.

- `source`:

  Source of the genome information: `"built-in"`, `"UCSC"`, or
  `"UCSC-chr.len"`.

- `cytoband`:

  Data frame with columns `chr`, `start`, `end`, `name`, and `gieStain`.

- `chromosomes`:

  Filtered chromosome vector in plotting order.

- `plot_type`:

  Default `circlize` plot type vector.

- `chr_lengths`:

  Named numeric vector of chromosome lengths.

## Details

This function is tolerant of different cytoband column names. It
standardizes cytoband data to `chr`, `start`, `end`, `name`, and
`gieStain` before plotting.

## Notes

For `ce10`, `ce11`, and `dm6`, synthetic cytobands are generated from
built-in chromosome lengths. These do not display real cytogenetic
bands; they provide chromosome axes and labels.

For UCSC-supported builds, the cytoband is obtained through
[`circlize::read.cytoband()`](https://rdrr.io/pkg/circlize/man/read.cytoband.html).
If only chromosome lengths are available, synthetic cytobands are
generated and the default plot type omits the ideogram.

## References

Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize implements and
enhances circular visualization in R. Bioinformatics.
2014;30(19):2811-2812.

## See also

[`nice_circos()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_circos.md),
[`is_circos_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/is_circos_genome.md)

## Examples

``` r
ce10 <- resolve_genome("ce10")
ce10
#> Genome: ce10 (built-in)
#> Chromosomes: chrI, chrII, chrIII, chrIV, chrV, chrX (6)
#> Plot type: axis, labels

resolve_genome("dm6", chromosome_index = c("chrX", "chr2L"))
#> Genome: dm6 (built-in)
#> Chromosomes: chrX, chr2L (2)
#> Plot type: axis, labels

if (FALSE) { # \dontrun{
hg38 <- resolve_genome("hg38", chromosome_index = "chr17")
print(hg38)
head(hg38$cytoband)
} # }
```
