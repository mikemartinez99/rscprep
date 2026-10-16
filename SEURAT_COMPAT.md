# Seurat 4 / Seurat 5 compatibility

rscprep now detects whether Seurat 4 or Seurat 5 is installed and
adjusts on its own. This note records what changed, which functions it
touches, and where the two versions actually behave differently.

Tested against Seurat 4.4.0 / SeuratObject 4.1.4 and Seurat 5.4.0 /
SeuratObject 5.3.0. The full test suite passes under both (676 tests
under 4, 682 under 5; the difference is tests that only apply to one
version and skip under the other).

------------------------------------------------------------------------

## Summary

| Function | Code changed? | Behaves differently between v4 and v5? |
|----|----|----|
| [`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md) | **Yes** | **Yes, the biggest difference.** Seurat 5’s [`merge()`](https://rdrr.io/r/base/merge.html) splits layers by sample. [`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md) now joins them by default. New argument `joinLayers`. |
| [`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md) | **Yes** | **Yes.** Raw-data detection was broken under Seurat 5 and is now fixed. |
| [`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md) | No | Only through `doCellCycle = TRUE`, which calls [`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md). |
| [`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md) | No | Returns `Assay5` objects under v5 and `Assay` under v4. Same counts, same metadata. |
| [`includeMetadata()`](https://mikemartinez99.github.io/rscprep/reference/includeMetadata.md) | No | No |
| [`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md), [`getOutlierSummary()`](https://mikemartinez99.github.io/rscprep/reference/getOutlierSummary.md) | No | No; they only read and write metadata. |
| [`rawMatrixFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/rawMatrixFrom10x.md) | No | No; it returns plain matrices, not Seurat objects. |
| [`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md) | Yes, for **DropletUtils**, not Seurat | See the separate DropletUtils section at the end. |
| [`theme_prepro()`](https://mikemartinez99.github.io/rscprep/reference/theme_prepro.md), `scale_*_prepro()`, [`prepro_pal()`](https://mikemartinez99.github.io/rscprep/reference/prepro_pal.md) | No | No |

------------------------------------------------------------------------

## How detection works

All version handling lives in one internal file, `R/seurat_compat.R`. No
exported function branches on the version itself; each one calls these
helpers.

| Helper | What it does |
|----|----|
| `isSeuratV5()` | `TRUE` when SeuratObject 5.0.0 or later is installed. It reads SeuratObject’s version, not Seurat’s, because the classes and accessors live in SeuratObject. |
| `seuratV5Fn(name)` | Fetches a SeuratObject-5-only function (`Layers`, `LayerData`, `JoinLayers`) at call time with [`getExportedValue()`](https://rdrr.io/r/base/ns-reflect.html). That way the package still installs and passes `R CMD check` against SeuratObject 4, where those functions don’t exist. |
| `getLayerValues(obj, layer, assay)` | Reads `"counts"`, `"data"` or `"scale.data"` the right way for each version: `GetAssayData(slot = )` on v4, `LayerData()` on v5. On v5 it also stitches split per-sample layers (`data.Sample_A`, `data.Sample_B`, …) back into one matrix, and returns `NULL` when the layer doesn’t exist. |
| `joinAssayLayers(obj)` | Seurat 5 only: joins the split layers of every `Assay5` assay. A no-op on Seurat 4, and for any classic `Assay` (which has no layers). |

`SeuratObject` was added to `Imports` in `DESCRIPTION`, since the
package now calls it directly. It was already installed in every case,
as a dependency of Seurat.

The object can be a different class from the installed version: an
object saved under Seurat 4 still loads under Seurat 5 as a classic
`Assay`, not an `Assay5`. So `joinAssayLayers()` checks each assay’s
class rather than relying on `isSeuratV5()` alone.

------------------------------------------------------------------------

## The three differences that matter

### 1. `merge()` returns a different shape (affects `mergeSeurat()`)

**This is the change most likely to affect downstream code.**

|  | Seurat 4 | Seurat 5 |
|----|----|----|
| [`merge()`](https://rdrr.io/r/base/merge.html) output | One `counts` matrix covering every cell | One layer per input: `counts.Sample_A`, `counts.Sample_B`, … |
| After `NormalizeData()` | One `data` matrix | `data.Sample_A`, `data.Sample_B`, … |

Code that expects one counts matrix (`GetAssayData()`,
`obj[["RNA"]]$counts`, anything that indexes the assay directly) gets a
different object depending on which version made it. Seurat 5’s own
integration workflow (`IntegrateLayers()`) needs the layers split; most
other code needs them joined.

**What changed in
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md):**

- **New argument `joinLayers = TRUE`.** Under Seurat 5 the merged
  object’s layers are joined, so you get one `counts` layer, the same
  shape Seurat 4 returns. Under Seurat 4 the argument does nothing.
- To keep the layers split for `IntegrateLayers()`, call
  `mergeSeurat(objList, joinLayers = FALSE)`.

**Behaviour change for existing Seurat 5 users:** before this change,
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
under Seurat 5 returned split layers. It now returns them joined. Any
Seurat 5 code that relied on the split, such as a call to
`IntegrateLayers()` straight after merging, needs `joinLayers = FALSE`.

**Bug fixed along the way:**
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
counted cells with `vapply(objList, ncol, integer(1))`. Under
SeuratObject 5, [`ncol()`](https://rdrr.io/r/base/nrow.html) returns a
double, so
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
errored on *every* Seurat 5 merge with
`values must be type 'integer', but FUN(X[[1]]) result is type 'double'`.
It now uses `numeric(1)`, which accepts both.

**New check, same message on both versions:**
[`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md)
now checks, before merging, that every object has an `Ident` metadata
column, and errors with
`The following element(s) of '<arg>' have no Ident column: <names>`.
Previously a missing `Ident` slipped through to the post-merge check,
where SeuratObject raised its own internal error. That error was worded
differently by v4 and v5 and named neither the argument nor the sample.

### 2. Unnormalised objects have no `data` layer (affects `runCellCycle()`)

|  | Seurat 4 | Seurat 5 |
|----|----|----|
| `data` right after `CreateSeuratObject()` | A copy of `counts` | **Doesn’t exist** until `NormalizeData()` runs |
| Asking for it | Returns the copy | `GetAssayData(layer = "data")` returns an empty 0 × 0 matrix with a warning |
| `GetAssayData(slot = "data")` | Works | **Error:** `slot` is defunct in SeuratObject 5 |

[`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md)
normalises raw data before scoring, because scoring assumes normalised
values. It decides whether data is raw with the internal
`looksUnnormalised()`, which checks whether the `data` values are all
whole numbers.

**What was broken under Seurat 5:** `looksUnnormalised()` called
`GetAssayData(slot = "data")`. Under SeuratObject 5 that errors. The
error was caught and turned into `NA` (“can’t tell”), so raw data was
**never detected**.
[`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md)
then skipped normalisation, scored raw counts, and usually failed inside
[`Seurat::CellCycleScoring()`](https://satijalab.org/seurat/reference/CellCycleScoring.html).

**What changed:** `looksUnnormalised()` now reads through
`getLayerValues()`:

- **v4:** reads the `data` slot, which holds the counts copy, and checks
  for whole numbers as before.
- **v5:** if there is no `data` layer at all, returns `TRUE` (definitely
  raw) without inspecting anything. If one exists, it checks for whole
  numbers as on v4.
- **v5, after a merge with split layers:** reads across all the `data.*`
  layers.

The same fix applies to `diagnoseQC(..., doCellCycle = TRUE)`, which
calls
[`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md).

[`Seurat::CellCycleScoring()`](https://satijalab.org/seurat/reference/CellCycleScoring.html)
itself was checked under Seurat 5 with split layers and works, so
[`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md)
does not need to join them first.

### 3. Addressing assay data: `slot` versus `layer`

SeuratObject 5 renamed `slot =` to `layer =` in `GetAssayData()` and
`SetAssayData()`, and has since made `slot` defunct. Any code calling
`GetAssayData(obj, slot = ...)` now errors under v5 rather than just
warning. In the package, the only place that did this was
`looksUnnormalised()` (above). It now goes through `getLayerValues()`.
The two remaining `slot =` calls in the codebase sit inside the Seurat 4
branches of `getLayerValues()` and the test helper `set_data_layer()`,
which is correct.

------------------------------------------------------------------------

## Functions checked and found version-neutral

These were not changed. Their output under v4 and v5 was compared
directly, by rendering the vignettes against both versions on real 10x
data during this work, and was identical apart from Seurat 5 printing an
extra `N layers present` line when an object is shown.

- **[`seuratFrom10x()`](https://mikemartinez99.github.io/rscprep/reference/seuratFrom10x.md):**
  `CreateSeuratObject()` builds a classic `Assay` under v4 and an
  `Assay5` under v5. Counts, `nCount_RNA`, `nFeature_RNA`, `orig.ident`
  and `Ident` are identical. The vignette used to say the function
  “forces classic (v4-style) Seurat assay objects internally”; it never
  did, and that sentence was removed in an earlier change. **The
  vignettes have since been reverted, so the sentence is back in
  `vignettes/seuratFrom10x.Rmd` and should be removed again.**
- **[`diagnoseQC()`](https://mikemartinez99.github.io/rscprep/reference/diagnoseQC.md):**
  [`Seurat::PercentageFeatureSet()`](https://satijalab.org/seurat/reference/PercentageFeatureSet.html)
  gives the same `pctMito` and `pctRibo` on joined and on split layers.
- **[`includeMetadata()`](https://mikemartinez99.github.io/rscprep/reference/includeMetadata.md),
  [`flagQC()`](https://mikemartinez99.github.io/rscprep/reference/flagQC.md),
  [`getOutlierSummary()`](https://mikemartinez99.github.io/rscprep/reference/getOutlierSummary.md):**
  these read and write `@meta.data` and `@misc` only, which did not
  change between versions.

------------------------------------------------------------------------

## Tests

| File | What it covers |
|----|----|
| `tests/testthat/test-seurat_compat.R` (new) | `isSeuratV5()`; `getLayerValues()` on counts, on a missing v5 `data` layer (`NULL`), and on split layers; `joinAssayLayers()` on both versions; [`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md) with `joinLayers = TRUE` and `FALSE`; [`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md) normalising raw data. Tests for Seurat 5-only behaviour skip under Seurat 4, and vice versa. |
| `tests/testthat/test-mergeSeurat.R` | The `Ident` test now expects the package’s own message rather than SeuratObject’s. |
| `tests/testthat/test-runCellCycle.R` | `GetAssayData(slot = )` calls replaced with `getLayerValues()`; `SetAssayData(slot = )` replaced with the new fixture `set_data_layer()`. |
| `tests/testthat/test-seuratFrom10x-manifest.R` | `vapply(..., integer(1))` changed to `numeric(1)`. |
| `tests/testthat/helper-fixtures.R` | `make_empty_seurat()` builds a classic `Assay` under either version (it writes `@counts`/`@data` directly, which an `Assay5` doesn’t have). New `set_data_layer()` writes normalised data under either version. |

To run the suite against Seurat 5 without touching your main library,
install Seurat 5 into a separate library and put it first on the path:

``` r

lib5 <- "~/R-seurat5"
dir.create(lib5)
install.packages(c("SeuratObject", "Seurat"), lib = lib5)

.libPaths(c(lib5, .libPaths()))
devtools::test()
```

------------------------------------------------------------------------

## Known gaps

- **Mixed lists.** A list that mixes classic-`Assay` objects (for
  example, one saved under Seurat 4) with `Assay5` objects made under
  Seurat 5 has not been tested through
  [`mergeSeurat()`](https://mikemartinez99.github.io/rscprep/reference/mergeSeurat.md).
- **Vignettes.** These were reverted to their state before the example
  data was added, so they predate this work:
  - `vignettes/seuratFrom10x.Rmd` has the incorrect “forces classic
    assay objects” sentence again.
  - `vignettes/plotBarcodeRanks.Rmd` recommends `df` and `fit.bounds`,
    which do nothing with DropletUtils 1.30 and later (see below).

  All vignettes are `eval = FALSE`, so none of this breaks the site
  build.
- **[`runCellCycle()`](https://mikemartinez99.github.io/rscprep/reference/runCellCycle.md)
  tests that use real ortholog conversion** skip when `orthogene` is not
  installed, so the non-human path has only been tested with the
  conversion mocked.

------------------------------------------------------------------------

## Separate: DropletUtils 1.30 (`plotBarcodeRanks()`)

This is not a Seurat change, but it came up in the same work and follows
the same pattern: the installed version is detected from what the
function returns.

From DropletUtils 1.30, `barcodeRanks()` finds the knee by tracing the
curve instead of fitting a smoothing spline. It no longer returns a
`fitted` column, and its `fit.bounds` and `df` arguments are ignored.
[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md)
read `fitted` unconditionally, so it failed with
`arguments imply differing number of rows: 20650, 0`.

**What changed in
[`plotBarcodeRanks()`](https://mikemartinez99.github.io/rscprep/reference/plotBarcodeRanks.md):**

- `data$fitted` is `NA` when DropletUtils returns no spline (1.30 and
  later), and holds the spline as before on older versions.
- The fitted line is drawn only when there is a spline to draw; the rest
  of the figure is unchanged.
- The help page has a new “DropletUtils versions” section:

| DropletUtils | How the knee is found | Arguments that tune it | Fitted line |
|----|----|----|----|
| before 1.28 | smoothing spline, then derivatives | `fit.bounds`, `df` | yes |
| 1.28 | spline; knee is the point farthest from the plateau-to-inflection line | `fit.bounds`, `df` | yes |
| 1.30 and later | curve tracing, no spline | `window`, `gradient.threshold` | no |

The same data can give a different knee, and sometimes a different
inflection, under different DropletUtils versions. Record
`packageVersion("DropletUtils")` alongside any cutoff you report.

Tests in `tests/testthat/test-plotBarcodeRanks.R` hand
`summariseBarcodeRanks()` both shapes of result directly, so both are
covered whichever DropletUtils is installed.
