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

## Item 2: Lightweight walk (no Cell class)

* Rewrote `scanline_burn.cpp` with lightweight walk — no Cell class allocation.
* `LightTraversal` struct tracks entry/exit coordinates and sides directly.
* `Box::crossing()` used directly for cell exit computation instead of
  `Cell::take()`.
* Coverage fractions computed via `left_hand_area()` (same algorithm as
  `Cell::covered_fraction`), with a shoelace shortcut for closed rings
  within a single cell.
* Fixed bug where zero-area boundary traversals (edges along cell walls)
  lost their winding information, causing missing interior cells for
  concave polygons and axis-aligned edges.
* The vendored `cell.h`, `cell.cpp`, `traversal.h`, `traversal.cpp` are
  no longer used by `burn_scanline()` — only by the original `burn_sparse()`.
* Validated against `burn_sparse()` — identical output (8/8 tests pass).

## Roadmap (remaining)

* Item 3: Analytical single-edge coverage (bypass left_hand_area for
  single-traversal cells — requires correct perimeter-distance logic to
  select the right CCW corner walk direction)
* Item 4: Benchmark perimeter-proportional scaling vs tiled dense approach
* Item 5: Multi-polygon shared-boundary handling
* Item 6: Edge cases (vertex on cell boundary, horizontal/vertical edges, slivers)
* Final: migrate to controlledburn as the production home
