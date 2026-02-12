# benchmark-scaling.R — Item 3: Benchmark perimeter-proportional scaling
#
# Compare burn_scanline (O(perimeter) memory, analytical coverage) against
# burn_sparse (O(bbox area) memory, dense matrix + flood fill).
#
# Key hypothesis: for a fixed polygon, doubling grid resolution 
# quadruples burn_sparse cost but only doubles burn_scanline cost
# (boundary cells scale with perimeter, interior runs are O(1) per row).

library(denseburn)
library(geos)

if (requireNamespace("bench", quietly = TRUE)) {
  has_bench <- TRUE
  library(bench)
} else {
  has_bench <- FALSE
  message("Install 'bench' for detailed benchmarks. Using system.time() fallback.")
}

# ---- Helper: time a single call ----
time_one <- function(expr, times = 3) {
  if (has_bench) {
    b <- bench::mark(expr, min_iterations = times, check = FALSE)
    as.numeric(b$median)  # seconds
  } else {
    timings <- replicate(times, system.time(expr)[["elapsed"]])
    median(timings)
  }
}

# ---- Test geometries ----

# 1. Simple square — compact, low perimeter:area
square_wkt <- "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))"

# 2. Thin horizontal sliver — high perimeter:area
sliver_wkt <- "POLYGON ((0.5 4.5, 9.5 4.5, 9.5 5.5, 0.5 5.5, 0.5 4.5))"

# 3. Star shape — complex boundary
star_coords <- function(n = 12, r_outer = 4.5, r_inner = 2.0, cx = 5, cy = 5) {
  angles <- seq(0, 2 * pi, length.out = 2 * n + 1)
  radii <- rep(c(r_outer, r_inner), length.out = 2 * n + 1)
  x <- cx + radii * cos(angles)
  y <- cy + radii * sin(angles)
  paste0("POLYGON ((", paste(sprintf("%.6f %.6f", x, y), collapse = ", "), "))")
}
star_wkt <- star_coords(12)

# 4. Donut — hole means more boundary per area
donut_wkt <- "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))"

# 5. Fractal-ish coastline — many vertices, high perimeter
jagged_coords <- function(n_teeth = 40, amplitude = 0.3) {
  # Sawtooth ring around a square
  base_x <- c(seq(1, 9, length.out = n_teeth),
               rep(9, n_teeth),
               seq(9, 1, length.out = n_teeth),
               rep(1, n_teeth))
  base_y <- c(rep(1, n_teeth),
               seq(1, 9, length.out = n_teeth),
               rep(9, n_teeth),
               seq(9, 1, length.out = n_teeth))
  # Add sawtooth perturbation perpendicular to edge direction
  n <- length(base_x)
  perturb <- rep(c(0, amplitude), length.out = n)
  # Determine perpendicular direction per segment
  px <- numeric(n); py <- numeric(n)
  side_n <- n_teeth
  # bottom edge: perturb in +y
  px[1:side_n] <- base_x[1:side_n]
  py[1:side_n] <- base_y[1:side_n] + perturb[1:side_n]
  # right edge: perturb in -x
  idx <- (side_n+1):(2*side_n)
  px[idx] <- base_x[idx] - perturb[idx]
  py[idx] <- base_y[idx]
  # top edge: perturb in -y
  idx <- (2*side_n+1):(3*side_n)
  px[idx] <- base_x[idx]
  py[idx] <- base_y[idx] - perturb[idx]
  # left edge: perturb in +x
  idx <- (3*side_n+1):(4*side_n)
  px[idx] <- base_x[idx] + perturb[idx]
  py[idx] <- base_y[idx]
  # Close the ring
  px <- c(px, px[1])
  py <- c(py, py[1])
  paste0("POLYGON ((", paste(sprintf("%.6f %.6f", px, py), collapse = ", "), "))")
}
jagged_wkt <- jagged_coords(40)

geoms <- list(
  square = as_geos_geometry(square_wkt),
  sliver = as_geos_geometry(sliver_wkt),
  star   = as_geos_geometry(star_wkt),
  donut  = as_geos_geometry(donut_wkt),
  jagged = as_geos_geometry(jagged_wkt)
)

# ---- Grid resolution sweep ----
#
# Scale grid from coarse to fine while holding extent fixed at (0, 10, 0, 10).
# This increases cell count as n^2 but perimeter cell count as n.

ext <- c(0, 10, 0, 10)
resolutions <- c(50, 100, 200, 400, 800, 1600)

cat("=== Scaling benchmark ===\n")
cat(sprintf("Resolutions: %s\n", paste(resolutions, collapse = ", ")))
cat(sprintf("Geometries: %s\n\n", paste(names(geoms), collapse = ", ")))

results <- list()

for (gname in names(geoms)) {
  g <- geoms[[gname]]
  cat(sprintf("--- %s ---\n", gname))
  cat(sprintf("%8s  %12s  %12s  %8s  %8s  %8s  %8s\n",
              "res", "scanline_ms", "sparse_ms", "ratio",
              "sl_runs", "sl_edges", "sp_runs"))

  for (res in resolutions) {
    dim <- c(res, res)

    t_scanline <- time_one(burn_scanline(g, extent = ext, dimension = dim)) * 1000
    t_sparse   <- time_one(burn_sparse(g, extent = ext, dimension = dim)) * 1000

    # Output sizes for structural comparison
    r_sl <- burn_scanline(g, extent = ext, dimension = dim)
    r_sp <- burn_sparse(g, extent = ext, dimension = dim)

    ratio <- t_sparse / t_scanline

    cat(sprintf("%8d  %12.2f  %12.2f  %8.2f  %8d  %8d  %8d\n",
                res, t_scanline, t_sparse, ratio,
                nrow(r_sl$runs), nrow(r_sl$edges),
                nrow(r_sp$runs)))

    results[[length(results) + 1]] <- data.frame(
      geometry = gname, resolution = res,
      scanline_ms = t_scanline, sparse_ms = t_sparse,
      ratio = ratio,
      sl_runs = nrow(r_sl$runs), sl_edges = nrow(r_sl$edges),
      sp_runs = nrow(r_sp$runs), sp_edges = nrow(r_sp$edges),
      stringsAsFactors = FALSE
    )
  }
  cat("\n")
}

results_df <- do.call(rbind, results)

# ---- Scaling analysis ----
cat("=== Scaling rates (log2 of time ratio between consecutive resolutions) ===\n")
cat("(O(n) ≈ 1.0, O(n²) ≈ 2.0, O(n² log n) ≈ 2.x)\n\n")

for (gname in unique(results_df$geometry)) {
  sub <- results_df[results_df$geometry == gname, ]
  if (nrow(sub) < 2) next

  sl_rates <- diff(log2(sub$scanline_ms))
  sp_rates <- diff(log2(sub$sparse_ms))
  res_labels <- paste0(sub$resolution[-1], "/", sub$resolution[-nrow(sub)])

  cat(sprintf("--- %s ---\n", gname))
  cat(sprintf("%12s  %10s  %10s\n", "step", "scanline", "sparse"))
  for (i in seq_along(sl_rates)) {
    cat(sprintf("%12s  %10.2f  %10.2f\n", res_labels[i], sl_rates[i], sp_rates[i]))
  }
  cat(sprintf("%12s  %10.2f  %10.2f\n", "median",
              median(sl_rates), median(sp_rates)))
  cat("\n")
}

# ---- Edge count scaling (should be ~linear in resolution for scanline) ----
cat("=== Edge count scaling (log2 ratio between consecutive resolutions) ===\n")
cat("(O(n) ≈ 1.0 — boundary cells grow linearly with resolution)\n\n")

for (gname in unique(results_df$geometry)) {
  sub <- results_df[results_df$geometry == gname, ]
  if (nrow(sub) < 2) next

  edge_rates <- diff(log2(pmax(sub$sl_edges, 1)))
  res_labels <- paste0(sub$resolution[-1], "/", sub$resolution[-nrow(sub)])

  cat(sprintf("%s: ", gname))
  for (i in seq_along(edge_rates)) {
    cat(sprintf("%.2f ", edge_rates[i]))
  }
  cat(sprintf(" (median: %.2f)\n", median(edge_rates)))
}

cat("\n=== Summary ===\n")
cat("If scanline scaling rate ≈ 1.0 and sparse ≈ 2.0, the O(perimeter)\n")
cat("hypothesis is confirmed: doubling resolution doubles scanline cost\n")
cat("but quadruples dense cost.\n")
