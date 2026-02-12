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

## Item 3: Analytical single-edge coverage

* For single-traversal boundary cells (~90% of all boundary cells), coverage
  fraction is now computed analytically: the traversal path + CW cell boundary
  corners from exit to entry form a closed polygon, area via shoelace formula.
  No `left_hand_area` chain-chasing, no heap allocation.
* Uses `perimeter_distance()` (same convention as `left_hand_area`) to
  determine which corners fall in the CW arc from exit to entry. This
  correctly handles all entry/exit side combinations including same-side cases.
* Multi-traversal cells (polygon vertex in cell, or multiple edges crossing)
  still fall back to `left_hand_area()`.
* Also added user-friendly defaults: `extent` derived from geometry bbox,
  `dimension` auto-fitted (256 cells on long axis), and `resolution` parameter
  as alternative to `dimension`. New dependency on `wk` (already a dependency
  of geos and sf).

## Item 4: Benchmark perimeter-proportional scaling

* `inst/docs-design/denseburn-refactor/benchmark-scaling.R`: scaling benchmark
  across grid resolutions (100–3200) for five geometry types.
* Confirmed O(perimeter) scaling for `burn_scanline`: median log2 rate ≈ 1.0
  (linear) across all geometry types. `burn_sparse` trends toward ≈ 2.0
  (quadratic) at higher resolutions.
* Speedup grows with resolution: star polygon 17× at 3200×3200, jagged
  coastline 9×, simple shapes 2–6×.
* Memory: 10M-cell grid requires 39 MB dense matrix vs 200–600 KB sparse
  output. Real-world test (CGAZ world polygons, 32K×16K = 500M cells):
  50 MB sparse vs ~2 GB dense.
* Edge count scaling confirms linear boundary growth: star/jagged median
  ≈ 1.0, axis-aligned shapes have zero edges (all boundaries on grid lines).

## Item 5: Multi-polygon shared-boundary validation

* `inst/docs-design/denseburn-refactor/validate-shared-boundary.R`:
  comprehensive test suite for coverage complementarity at shared boundaries.
* Verified that adjacent polygons produce coverage fractions summing to
  exactly 1.0 at every shared cell (18/18 tests pass).
* Test cases: vertical/horizontal shared edges (on-grid and off-grid),
  four quadrants meeting at a point, diagonal shared edges (two triangles
  forming a square), irregular tilings (L-shape + complement), donut with
  inner fill, high-res diagonal (100×100), and NC counties (real data).
* No code changes needed: `Box::crossing()` is deterministic for coincident
  edges, so independent per-polygon walks produce complementary fractions.

## Roadmap (remaining)

* Item 6: Edge cases (vertex on cell boundary, horizontal/vertical edges, slivers)
* Final: migrate to controlledburn as the production home
