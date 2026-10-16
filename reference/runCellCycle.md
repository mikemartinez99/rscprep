# Score Cells by Cell Cycle Phase

Assigns each cell an S-phase score, a G2/M score and a phase call, using
the Tirosh et al. marker sets that ship with Seurat and
[`Seurat::CellCycleScoring()`](https://satijalab.org/seurat/reference/CellCycleScoring.html).
The marker lists are human symbols; for any other organism they are
converted to that species' orthologs first, so the same call works
across species without the caller assembling gene lists.

## Usage

``` r
runCellCycle(
  obj,
  organism = c("hsa", "mmu", "rno", "dre"),
  method = c("gprofiler", "homologene", "babelgene")
)
```

## Arguments

- obj:

  A single
  [Seurat](https://satijalab.org/seurat/reference/reexports.html)
  object, normalised.

- organism:

  Three-letter code selecting the marker set: `"hsa"` (human, the
  default) uses
  [Seurat::cc.genes](https://satijalab.org/seurat/reference/cc.genes.html)
  directly; `"mmu"` (mouse), `"rno"` (rat) and `"dre"` (zebrafish)
  convert it to orthologs first.

- method:

  Ortholog mapping backend passed to
  [`orthogene::convert_orthologs()`](https://rdrr.io/pkg/orthogene/man/convert_orthologs.html).
  `"gprofiler"` (the default) is the most complete but needs network
  access; `"homologene"` and `"babelgene"` resolve locally. Ignored when
  `organism = "hsa"`.

## Value

The input object with `S.Score`, `G2M.Score` and `Phase` added to its
metadata. Cell count and identities are unchanged. The data layer is
unchanged too, unless the object arrived unnormalised – see
*Normalisation*.

## Details

Cell cycle is a QC concern because it is one of the largest sources of
variation in a proliferating population, and it is *structured*
variation – cycling cells cluster with cycling cells regardless of type.
The scores are added here so that they can be looked at, and if
necessary regressed out later; nothing is regressed, scaled or clustered
by this function.

## What gets added

[`Seurat::CellCycleScoring()`](https://satijalab.org/seurat/reference/CellCycleScoring.html)
writes three metadata columns:

- `S.Score`, `G2M.Score`:

  Mean expression of each marker set, centred against a background of
  control features drawn from comparable expression bins. Relative, not
  absolute: a score near zero means "no more than the background", not
  "no expression".

- `Phase`:

  `"S"`, `"G2M"` or `"G1"`, from whichever score is higher, with `"G1"`
  when both are negative. `"G1"` therefore also absorbs G0 – the two are
  not separable from these markers.

`Idents()` is deliberately left alone (`set.ident = FALSE`), so the
sample identity
[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md)
set survives.

## Normalisation

Scoring bins genes by average expression, which assumes the data have
been normalised. An object that still looks like raw counts is
normalised here with
[`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html)
at its defaults, after a warning – and the returned object carries that
normalised data, since scoring against one version of the data and
handing back another would be worse than either.

Normalise beforehand if you want anything other than the defaults; the
check is a heuristic (log-normalised values are essentially never whole
numbers), so an already-normalised object is left exactly as it was.

The installed Seurat version is detected automatically. Seurat 4 copies
the counts into the data slot when an object is created, so the
heuristic reads that copy. A Seurat 5 object has no data layer until it
is normalised, so its absence alone means raw. Seurat 5 objects whose
layers are still split by sample after a merge are read across every
layer.

## Non-human organisms

Markers are converted with
[`orthogene::convert_orthologs()`](https://rdrr.io/pkg/orthogene/man/convert_orthologs.html),
which is only needed – and only loaded – when `organism` is not `"hsa"`.
Genes without a one-to-one ortholog are dropped, so the converted sets
are smaller than the human ones; how many survived is reported.

The default conversion `method`, `"gprofiler"`, queries a web service
and so needs network access. On a machine without it, pass
`method = "homologene"`, which resolves locally from a bundled mapping
at the cost of some coverage.

## See also

[`Seurat::CellCycleScoring()`](https://satijalab.org/seurat/reference/CellCycleScoring.html)
for the scoring itself,
[`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md)
for the other per-cell QC metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat <- mergeSeurat(seuratFrom10x(path = "/path/to/Pool1/outs", multi = TRUE))
seurat <- Seurat::NormalizeData(seurat)

seurat <- runCellCycle(seurat, organism = "mmu")

table(seurat$Phase)

# is one sample cycling harder than the others?
table(seurat$Ident, seurat$Phase)

# on a cluster with no outbound network
seurat <- runCellCycle(seurat, organism = "mmu", method = "homologene")

# the scores are ordinary metadata, so they plot like any other metric
plotScatter(seurat, "S.Score", "G2M.Score", group_by = "Phase")
} # }
```
