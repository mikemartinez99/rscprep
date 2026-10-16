# Merge a Named List of Seurat Objects into One

Combines a named list of per-sample
[Seurat](https://satijalab.org/seurat/reference/reexports.html) objects
– typically the output of
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
– into a single object via
[`merge()`](https://rdrr.io/r/base/merge.html), prefixing every cell
barcode with its sample name so barcodes stay unique and traceable back
to their sample of origin.

## Usage

``` r
mergeSeurat(objList)
```

## Arguments

- objList:

  Named list of
  [Seurat](https://satijalab.org/seurat/reference/reexports.html)
  objects, one per sample, as returned by
  [`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md).
  Names are used as cell barcode prefixes and must be unique.

## Value

A single merged
[Seurat](https://satijalab.org/seurat/reference/reexports.html) object,
returned **invisibly** (the object is printed as a side effect). Cell
barcodes are prefixed with their sample name; the `Ident` column is
carried through from the inputs. If `objList` has exactly one element,
that element is returned unmerged and unprefixed, with a message.

## Details

The list names are load-bearing: they are passed to
[`merge()`](https://rdrr.io/r/base/merge.html) as `add.cell.id`, so a
cell `AAACCCAAGA-1` in list entry `Sample_1` becomes
`Sample_1_AAACCCAAGA-1` in the merged object. Names must therefore be
present, non-blank and unique.

The whole list is validated up front, so a malformed batch fails with a
specific message instead of partway through the merge. In order,
`objList` must:

- be non-empty (a single-element list short-circuits – see *Value*);

- contain only structurally valid Seurat objects;

- be fully named, with no `NA`, blank or duplicated names;

- contain no zero-cell objects;

- agree on
  [`Seurat::DefaultAssay()`](https://satijalab.org/seurat/reference/reexports.html)
  across every element.

Two further checks run *after* the merge, guarding against a merge that
silently loses data: every sample must still be represented in the
merged `Ident` metadata column, and the merged cell count must equal the
sum of the inputs'. The first of these relies on each input object
carrying an `Ident` column, which
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
sets to the sample name.

## See also

[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md),
which produces the named list this expects.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_list <- seuratFrom10x(
  path        = "/path/to/Pool1/outs",
  multi = TRUE
)

merged <- mergeSeurat(seurat_list)

# barcodes are namespaced by sample
head(colnames(merged))

# and Ident still identifies the sample of origin
table(merged$Ident)
} # }
```
