#' Scanline polygon rasterization with exact coverage fractions
#'
#' Experimental scanline sweep implementation of sparse polygon rasterization.
#' Produces identical output to [burn_sparse()] but uses a winding-number
#' scanline sweep instead of flood fill for interior classification.
#' Memory usage is O(perimeter) rather than O(bounding-box area).
#'
#' @inheritParams burn_sparse
#'
#' @return A list with class `"denseburn"` containing:
#'   \describe{
#'     \item{`runs`}{data.frame with columns `row`, `col_start`, `col_end`, `id`}
#'     \item{`edges`}{data.frame with columns `row`, `col`, `weight`, `id`}
#'     \item{`extent`}{the raster extent as supplied}
#'     \item{`dimension`}{the grid dimensions as supplied}
#'   }
#'
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   library(geos)
#'   poly <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
#'   result <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
#'   result$runs
#'   result$edges
#' }
burn_scanline <- function(x, extent, dimension) {
  wkb <- as_wkb_list(x)
  extent <- as.double(extent)
  dimension <- as.integer(dimension)

  stopifnot(
    length(extent) == 4,
    length(dimension) == 2,
    dimension[1] > 0,
    dimension[2] > 0,
    extent[2] > extent[1],
    extent[4] > extent[3]
  )

  ncol_full <- dimension[1]
  nrow_full <- dimension[2]

  result <- cpp_scanline_burn(
    wkb,
    extent[1], extent[3], extent[2], extent[4],
    ncol_full, nrow_full
  )

  result$extent <- extent
  result$dimension <- dimension
  class(result) <- "denseburn"
  result
}
