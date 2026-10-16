# Add Per-Cell QC Metrics to a Seurat Object

Calculates the three per-cell metrics that single-cell filtering
decisions are usually made from – mitochondrial content, ribosomal
content and library complexity – and writes them into the object's
metadata as `pctMito`, `pctRibo` and `log10GenesPerUMI`. Nothing is
filtered and no cells are dropped; the metrics are added so that
thresholds can be chosen by looking at the distributions first.

## Usage

``` r
diagnoseQC(obj, organism = c("hsa", "mmu", "rno", "dre"), doCellCycle = FALSE)
```

## Arguments

- obj:

  A single
  [Seurat](https://satijalab.org/seurat/reference/reexports.html)
  object, typically the merged output of
  [`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md).
  A named *list* of objects is not accepted – merge it first, or
  [`lapply()`](https://rdrr.io/r/base/lapply.html) over it – so that the
  metrics are always calculated on the same object the QC plots are
  drawn from.

- organism:

  Three-letter code selecting the gene patterns: `"hsa"` (human, the
  default), `"mmu"` (mouse), `"rno"` (rat) or `"dre"` (zebrafish). See
  *Gene patterns by organism*.

- doCellCycle:

  Also score cell cycle phase, by handing the annotated object to
  [`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md)
  with the same `organism`. `FALSE` by default: unlike the metrics
  above, scoring assumes normalised data, and will normalise the object
  itself if it arrives raw. For a non-human organism it also needs
  orthogene. Call
  [`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md)
  directly if you need its `method` argument.

## Value

The input object with three numeric metadata columns added: `pctMito`,
`pctRibo` and `log10GenesPerUMI`. Existing columns of the same name are
overwritten, so the function is safe to re-run – with a different
`organism`, for instance. Every other slot is untouched, and the cell
count is unchanged.

With `doCellCycle = TRUE`, `S.Score`, `G2M.Score` and `Phase` are added
as well – see
[`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md)
for what they mean.

## Details

The two percentage metrics are computed by
[`Seurat::PercentageFeatureSet()`](https://satijalab.org/seurat/reference/PercentageFeatureSet.html)
from organism-specific gene *symbol* patterns, so the object's rownames
must be symbols. A matrix still keyed by Ensembl IDs will match nothing
and produce columns of zeroes; that case is warned about rather than
silently returned, and the warning names the pattern that failed so it
is obvious what went wrong.

## The metrics

- `pctMito`:

  Percentage of a cell's UMIs coming from mitochondrial genes. High
  values mark cells that lost cytoplasmic RNA through a ruptured
  membrane – dying or stressed cells – and are the usual first filter.
  What counts as high is tissue-dependent: 5-10% is a common starting
  point, but cardiac, muscle and kidney samples are legitimately
  mitochondria-rich and a fixed cutoff will throw away real cells.

- `pctRibo`:

  Percentage coming from ribosomal protein genes. Less often used as a
  hard filter and more often as a covariate: it varies strongly by cell
  type, so a sharp shift between samples that should be comparable
  points at a technical difference rather than a biological one.

- `log10GenesPerUMI`:

  Library complexity – `log10(nFeature) / log10(nCount)`. Falls between
  0 and 1 for real cells, and low values mean a library dominated by few
  transcripts: ambient RNA, degraded cells, or red blood cells, which
  are genuinely low-complexity. Around 0.8 is a conventional floor.

The count columns the complexity metric is derived from
(`nCount_<assay>` / `nFeature_<assay>`) are read from the object's
default assay, so an object whose default assay is not `RNA` is handled
without any extra argument. The three metrics themselves are always
named the same way, because that is what the plotting functions look
for.

## Gene patterns by organism

- `"hsa"` (human):

  mito `^MT-`, ribo `^RP[SL]`

- `"mmu"` (mouse):

  mito `^mt`, ribo `^Rp[sl]`

- `"rno"` (rat):

  mito `^Mt-`, ribo `^Rp[sl]`

- `"dre"` (zebrafish):

  mito `^mt-`, ribo `^rp[sl]`

These differ only in the capitalisation convention each species' gene
symbols follow, so an organism not listed here can borrow the code of
one whose convention it shares.

## See also

[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md),
which produces the object this expects;
[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md)
for the other side of QC – cell calling, from the raw matrix;
[`Seurat::PercentageFeatureSet()`](https://satijalab.org/seurat/reference/PercentageFeatureSet.html)
for the percentage calculation itself.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_list <- seuratFrom10x(path = "/path/to/Pool1/outs", multi = TRUE)
seurat      <- mergeSeurat(seurat_list)

seurat <- diagnoseQC(seurat, organism = "mmu")
head(seurat@meta.data)

# or metric-by-sample, without merging
seurat_list <- lapply(seurat_list, diagnoseQC, organism = "mmu")
} # }
```
