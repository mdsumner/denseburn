## compare-controlledburn.R — cross-package comparison
##
## Compares denseburn::burn_scanline() against controlledburn::burn_polygon()
## to verify agreement on interior cell classification.
##
## controlledburn does approximate cell-centre containment (binary in/out).
## denseburn does exact coverage fractions with boundary cells.
##
## Expected: interior runs (coverage = 1.0) should largely agree.
## Boundary cells differ: controlledburn classifies as in or out based on
## cell centre; denseburn gives the exact fraction.
##
## Requires: controlledburn, denseburn, sf (for NC test data)

if (!requireNamespace("controlledburn", quietly = TRUE)) {
  stop("Install controlledburn: remotes::install_github('hypertidy/controlledburn')")
}

library(denseburn)
library(controlledburn)

## ---- Helper: convert controlledburn output to dense matrix ----
cb_to_matrix <- function(cb_result, dm) {
  # cb_result is a list of 4-vectors: (start, end, row, poly_id), 0-based
  index <- matrix(unlist(cb_result, use.names = FALSE), ncol = 4L, byrow = TRUE) + 1L
  mat <- matrix(0, nrow = dm[2], ncol = dm[1])
  for (i in seq_len(nrow(index))) {
    mat[index[i, 3], index[i, 1]:index[i, 2]] <- 1
  }
  mat
}

## ---- Helper: compare ----
compare_cb_db <- function(sf_obj, ext, dm, label = "") {
  # controlledburn
  cb <- controlledburn:::burn_polygon(sf_obj, extent = ext, dimension = dm)
  mat_cb <- cb_to_matrix(cb, dm)

  # denseburn (using sf input)
  geom <- sf::st_geometry(sf_obj)
  db <- burn_scanline(geom, extent = ext, dimension = dm)
  mat_db <- materialise_chunk(db)

  # Interior cells in denseburn (coverage >= 0.999)
  mat_db_interior <- ifelse(mat_db >= 0.999, 1, 0)

  # Agreement on interior classification
  agree_interior <- sum(mat_cb == 1 & mat_db_interior == 1)
  cb_only <- sum(mat_cb == 1 & mat_db_interior == 0)
  db_only <- sum(mat_cb == 0 & mat_db_interior == 1)

  # Boundary cells in denseburn (0 < coverage < 1)
  n_boundary <- sum(mat_db > 0.001 & mat_db < 0.999)

  # Cells where controlledburn says "in" but denseburn says partial
  cb_in_db_partial <- sum(mat_cb == 1 & mat_db > 0.001 & mat_db < 0.999)
  # Cells where controlledburn says "out" but denseburn has partial coverage
  cb_out_db_partial <- sum(mat_cb == 0 & mat_db > 0.001 & mat_db < 0.999)

  cat(sprintf("--- %s (%dx%d) ---\n", label, dm[1], dm[2]))
  cat(sprintf("  Interior agreement:     %d cells\n", agree_interior))
  cat(sprintf("  cb=in, db=interior:     %d (cb extra interior)\n", cb_only))
  cat(sprintf("  cb=out, db=interior:    %d (db extra interior)\n", db_only))
  cat(sprintf("  denseburn boundary:     %d cells with 0 < coverage < 1\n", n_boundary))
  cat(sprintf("  cb=in, db=partial:      %d (cb calls partial cells interior)\n", cb_in_db_partial))
  cat(sprintf("  cb=out, db=partial:     %d (cb calls partial cells exterior)\n", cb_out_db_partial))

  # The total discrepancy should be small and concentrated at boundaries
  total_discrepancy <- cb_only + db_only
  cat(sprintf("  Total interior disagreement: %d cells\n\n", total_discrepancy))

  invisible(list(
    agree = agree_interior, cb_only = cb_only, db_only = db_only,
    boundary = n_boundary, cb_in_partial = cb_in_db_partial,
    cb_out_partial = cb_out_db_partial
  ))
}


## ---- Test data ----
cat("=== controlledburn vs denseburn comparison ===\n\n")

# NC counties — real-world polygons
if (requireNamespace("sf", quietly = TRUE)) {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  ext <- as.numeric(sf::st_bbox(nc))[c(1, 3, 2, 4)]

  compare_cb_db(nc[1:5, ], ext, c(200L, 80L), "NC 5 counties, low res")
  compare_cb_db(nc[1:5, ], ext, c(500L, 200L), "NC 5 counties, med res")
  compare_cb_db(nc[1:5, ], ext, c(1000L, 400L), "NC 5 counties, high res")

  compare_cb_db(nc, ext, c(500L, 200L), "NC all counties, med res")

  # Timing comparison
  cat("=== Timing comparison ===\n\n")
  geom <- sf::st_geometry(nc)
  dm <- c(2000L, 800L)

  cat(sprintf("Grid: %dx%d\n", dm[1], dm[2]))
  cat("controlledburn: ")
  print(system.time(cb <- controlledburn:::burn_polygon(nc, extent = ext, dimension = dm)))
  cat("denseburn (scanline): ")
  print(system.time(db <- burn_scanline(geom, extent = ext, dimension = dm)))
  cat("denseburn (sparse): ")
  print(system.time(db2 <- burn_sparse(geom, extent = ext, dimension = dm)))

  cat(sprintf("\nOutput sizes:\n"))
  cat(sprintf("  controlledburn: %s\n", format(object.size(cb), units = "auto")))
  cat(sprintf("  denseburn:      %s\n", format(object.size(db), units = "auto")))
}

cat("\n=== Notes ===\n")
cat("Interior disagreements are expected at polygon boundaries where\n")
cat("controlledburn uses cell-centre containment (binary) and denseburn\n")
cat("computes exact coverage fractions. The number of disagreements should\n")
cat("be O(perimeter) — proportional to the number of boundary cells.\n")
