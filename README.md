# denseburn

Sparse polygon rasterization with exact coverage fractions.

**denseburn** is the development repo for refactoring the dense-matrix
approach in [gridburn](https://github.com/hypertidy/gridburn) into a
pure scanline sweep that is O(perimeter) in memory. The final production
home will be [controlledburn](https://github.com/hypertidy/controlledburn).

## Status

| Item | Description | Status |
|------|-------------|--------|
| 1 | Scanline sweep prototype (Cell class + winding count) | ✓ validated |
| 2 | Lightweight walk (no Cell class, Box::crossing direct) | ✓ validated |
| 3 | Analytical single-edge coverage (perimeter_distance) | ✓ validated |
| 4 | Benchmark perimeter-proportional scaling | ✓ confirmed O(n) |
| 5 | Multi-polygon shared-boundary handling | ✓ complementary |
| 6 | Edge cases (vertex on boundary, horizontal edges, slivers) | planned |

## Installation

```r
# install.packages("pak")
pak::pak("hypertidy/denseburn")
```

## Usage

Two algorithms, identical output format:

```r
library(denseburn)
library(geos)

poly <- as_geos_geometry(
  "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))"
)
ext <- c(0, 10, 0, 10)
dim <- c(20L, 20L)

# Original: dense matrix + flood fill (from gridburn)
r1 <- burn_sparse(poly, extent = ext, dimension = dim)

# New: scanline winding sweep, O(perimeter) memory
r2 <- burn_scanline(poly, extent = ext, dimension = dim)

# Identical output
all.equal(materialise_chunk(r1), materialise_chunk(r2))
```

## Output format

Both functions return a list with two data.frames:

- **`runs`**: `(row, col_start, col_end, id)` — run-length encoded interior
  cells (coverage ≈ 1.0)
- **`edges`**: `(row, col, weight, id)` — boundary cells with exact partial
  coverage (0 < weight < 1)

## Lineage

```
fasterize (Wylie et al. scanline, approximate, dense output)
    └── controlledburn (approximate, sparse output, no pixel materialisation)
            └── gridburn (exact via vendored exactextract, sparse output, dense intermediate)
                    └── denseburn (exact, sparse output, no dense intermediate)
                            └── controlledburn (final home, both approximate and exact modes)
```

## License

Apache License 2.0
