library(denseburn)
library(geos)

poly <- as_geos_geometry(
  silicate::inlandwaters
)

ext <- as.numeric(wk::wk_bbox(poly))[c(1, 3, 2, 4)]
dim <- c(100L, 80L)

# Original: dense matrix + flood fill (from gridburn)
system.time(r1 <- burn_sparse(poly, extent = ext, dimension = dim))

# New: scanline winding sweep, O(perimeter) memory
system.time(r2 <- burn_scanline(poly, extent = ext, dimension = dim))


ximage::ximage(denseburn::mat)


system.time(x <- geos::as_geos_geometry(wk::wkb(vapour::vapour_read_geometry(sds::CGAZ()))))
system.time(r1 <- burn_scanline(x))

system.time(ximage::ximage(materialise_chunk(r1), r1$extent, col = hcl.colors(25)))



user  system elapsed
1.191   0.001   1.191

user  system elapsed
0.062   0.000   0.062
