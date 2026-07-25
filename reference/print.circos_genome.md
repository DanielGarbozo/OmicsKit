# Print a circos genome descriptor

Print a circos genome descriptor

## Usage

``` r
# S3 method for class 'circos_genome'
print(x, ...)
```

## Arguments

- x:

  A `"circos_genome"` object returned by
  [`resolve_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/resolve_genome.md).

- ...:

  Additional arguments, currently ignored.

## Value

Invisibly returns `x`.

## See also

[`resolve_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/resolve_genome.md)

## Examples

``` r
ce10 <- resolve_genome("ce10")
print(ce10)
#> Genome: ce10 (built-in)
#> Chromosomes: chrI, chrII, chrIII, chrIV, chrV, chrX (6)
#> Plot type: axis, labels
```
