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


ximage::ximage(materialise_chunk(r1), r1$extent)
ximage::ximage(materialise_chunk(r1), r2$extent)


system.time(x <- geos::as_geos_geometry(wk::wkb(vapour::vapour_read_geometry(sds::CGAZ()))))
system.time(r1 <- burn_scanline(x))
system.time(ximage::ximage(materialise_chunk(r1), r1$extent, col = hcl.colors(25)))


system.time(r1 <- burn_scanline(x, dimension = c(8192L, 4096L)))

system.time(r1 <- burn_scanline(x, dimension = c(8192L, 4096L) * 4))

## at item 2 we saw
#system.time(r1 <- burn_scanline(x, dimension = c(8192, 4096)))
#user  system elapsed
#1.890   0.000   1.889
#> pryr::object_size(r1)
#10.69 MB
#> system.time(r1 <- burn_scanline(x, dimension = c(8192, 4096) * 4))
#user  system elapsed
#3.703   0.000   3.701
#> pryr::object_size(r1)
#50.33 MB



