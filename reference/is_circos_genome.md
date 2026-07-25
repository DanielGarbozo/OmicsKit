# Test whether an object is a circos genome descriptor

Test whether an object is a circos genome descriptor

## Usage

``` r
is_circos_genome(x)
```

## Arguments

- x:

  An R object.

## Value

Logical scalar. `TRUE` if `x` inherits from class `"circos_genome"`,
otherwise `FALSE`.

## See also

[`resolve_genome()`](https://danielgarbozo.github.io/OmicsKit/reference/resolve_genome.md),
[`nice_circos()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_circos.md)

## Examples

``` r
ce10 <- resolve_genome("ce10")
is_circos_genome(ce10)
#> [1] TRUE
is_circos_genome(data.frame())
#> [1] FALSE
```
