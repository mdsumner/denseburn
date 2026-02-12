# denseburn 0.1.0

Forked from gridburn as the development home for the scanline refactor.

## Initial state (from gridburn)

* `burn_sparse()` — original dense-matrix approach using vendored exactextract
  algorithm (flood fill for interior classification, dense `Matrix<float>` for
  coverage fractions). Works correctly, validated against exactextractr. Tiled
  processing for large grids.

## Item 1: Scanline sweep prototype

* `burn_scanline()` — new function implementing the winding-number scanline
  sweep from the denseburn-refactor design. Reuses exactextract's `Cell` class
  for exact boundary coverage fractions, replaces flood fill with per-row
  winding count. Memory is O(perimeter) not O(bounding-box area).
* Validated cell-for-cell against `burn_sparse()` across simple rectangles,
  triangles, diamonds, L-shapes, polygons with holes, multi-polygons, and NC
  counties (8/8 tests pass).
* Does not yet use the vendored `floodfill.h`, `raster_cell_intersection.h`,
  `matrix.h`, or `raster.h` — these are only needed by the original
  `burn_sparse()` path.

## Roadmap (in this repo)

* Item 2: Analytical single-edge coverage (bypass Cell class for common case)
* Item 3: Benchmark perimeter-proportional scaling vs tiled dense approach
* Item 4: Multi-polygon shared-boundary handling
* Item 5: Edge cases (vertex on cell boundary, horizontal/vertical edges, slivers)
* Final: migrate to controlledburn as the production home
