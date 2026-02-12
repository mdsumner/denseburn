// scanline_burn.cpp — scanline polygon rasterization with exact coverage fractions
//
// This is the "controlledburn" algorithm from the denseburn-refactor design:
// reuse exactextract's Cell class for exact boundary cell coverage fractions,
// but replace the flood fill with a scanline winding-number sweep for interior
// classification. Memory usage is O(perimeter), not O(bounding-box area).
//
// The output is identical to gridburn's cpp_burn_sparse: two tables of
// (runs, edges) in sparse format.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include <cpp11.hpp>
#include <cpp11/list.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/raws.hpp>
#include <cpp11/strings.hpp>

#include <cstdarg>
#include <map>
#include <algorithm>
#include <cmath>

#include "libgeos.h"
#include "dense_to_sparse.h"

#include "exactextract/cell.h"
#include "exactextract/grid.h"
#include "exactextract/box.h"
#include "exactextract/geos_utils.h"
#include "exactextract/coordinate.h"
#include "exactextract/side.h"

// ---- GEOS context management (same as gridburn.cpp) ----

static void sl_geos_notice_handler(const char* /*fmt*/, ...) {}
static void sl_geos_error_handler(const char* fmt, ...) {
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    cpp11::stop("GEOS error: %s", buf);
}

class SLGEOSContextGuard {
public:
    SLGEOSContextGuard() {
        ctx_ = GEOS_init_r();
        if (!ctx_) cpp11::stop("Failed to create GEOS context");
        GEOSContext_setNoticeHandler_r(ctx_, sl_geos_notice_handler);
        GEOSContext_setErrorHandler_r(ctx_, sl_geos_error_handler);
    }
    ~SLGEOSContextGuard() { if (ctx_) GEOS_finish_r(ctx_); }
    GEOSContextHandle_t get() const { return ctx_; }
    SLGEOSContextGuard(const SLGEOSContextGuard&) = delete;
    SLGEOSContextGuard& operator=(const SLGEOSContextGuard&) = delete;
private:
    GEOSContextHandle_t ctx_;
};

class SLGEOSGeomGuard {
public:
    SLGEOSGeomGuard(GEOSContextHandle_t ctx, GEOSGeometry* g) : ctx_(ctx), g_(g) {}
    ~SLGEOSGeomGuard() { if (g_) GEOSGeom_destroy_r(ctx_, g_); }
    GEOSGeometry* get() const { return g_; }
    bool valid() const { return g_ != nullptr; }
    SLGEOSGeomGuard(const SLGEOSGeomGuard&) = delete;
    SLGEOSGeomGuard& operator=(const SLGEOSGeomGuard&) = delete;
private:
    GEOSContextHandle_t ctx_;
    GEOSGeometry* g_;
};

// ---- Per-cell boundary data for the scanline sweep ----

struct BoundaryCellRecord {
    int col;               // 0-based column in full grid
    float coverage;        // accumulated coverage fraction (signed: +exterior, -hole)
    int winding_delta;     // accumulated winding contribution
};

// ---- Scanline algorithm ----

// Walk a polygon ring through grid cells, collecting boundary cell data.
//
// This reuses the same Cell-based walk as exactextract's process_line, but
// instead of writing to a dense matrix, we record boundary cell coverage
// and winding crossings into per-row vectors.
//
// coords: ring vertices (will be normalised to CCW internally)
// is_ccw: whether the original ring is counter-clockwise
// is_exterior: true for exterior ring, false for hole
// grid: the infinite_extent grid for cell lookups
// full_grid: the bounded_extent grid for row/col mapping
// dy: cell height
// row_data: output — per-row vectors of boundary cell records (indexed by subgrid row)
// sub_rows, sub_cols: dimensions of the subgrid (geometry bbox within full grid)
// row_off, col_off: offset of subgrid in full grid
static void walk_ring(
    std::vector<exactextract::Coordinate> coords,
    bool is_ccw,
    bool is_exterior,
    const exactextract::Grid<exactextract::infinite_extent>& grid,
    double dy,
    std::vector<std::vector<BoundaryCellRecord>>& row_data,
    size_t sub_rows,
    size_t sub_cols,
    size_t row_off,
    size_t col_off
) {
    using namespace exactextract;

    if (coords.size() < 4) return; // degenerate ring

    // Normalise to CCW for correct Cell coverage fractions
    if (!is_ccw) {
        std::reverse(coords.begin(), coords.end());
    }

    // The sign factor: exterior rings contribute positively, holes negatively
    // (same as exactextract's add_ring_results)
    float coverage_factor = is_exterior ? 1.0f : -1.0f;
    // Winding factor: for CCW exterior rings, winding is standard (+1 up, -1 down).
    // For CW holes that we reversed to CCW, we negate the winding.
    int winding_factor = is_exterior ? 1 : -1;

    // Lazy Cell allocation: map from (row, col) in subgrid coords to Cell
    // (subgrid row/col are 1-based in the infinite_extent grid, but we use
    // the infinite_extent grid's row/col directly)
    std::map<std::pair<size_t, size_t>, std::unique_ptr<Cell>> cells;

    auto get_cell = [&](size_t r, size_t c) -> Cell* {
        auto key = std::make_pair(r, c);
        auto it = cells.find(key);
        if (it == cells.end()) {
            auto box = grid_cell(grid, r, c);
            auto result = cells.emplace(key, std::make_unique<Cell>(box));
            return result.first->second.get();
        }
        return it->second.get();
    };

    // ---- Walk boundary through cells ----
    // This replicates the core loop from raster_cell_intersection.cpp process_line

    size_t pos = 0;
    size_t row = grid.get_row(coords.front().y);
    size_t col = grid.get_column(coords.front().x);
    const Coordinate* last_exit = nullptr;

    // Track traversal entry/exit for winding calculation.
    // We record one entry per cell visit (one traversal).
    struct TraversalWinding {
        size_t row;
        size_t col;
        double entry_y;
        double exit_y;
    };
    std::vector<TraversalWinding> traversal_windings;

    while (pos < coords.size()) {
        Cell& cell = *get_cell(row, col);

        while (pos < coords.size()) {
            const Coordinate* next_coord = last_exit ? last_exit : &coords[pos];
            const Coordinate* prev_coord = pos > 0 ? &coords[pos - 1] : nullptr;

            cell.take(*next_coord, prev_coord);

            if (cell.last_traversal().exited()) {
                const Coordinate& exc = cell.last_traversal().exit_coordinate();
                if (exc != *next_coord) {
                    last_exit = &exc;
                }
                break;
            } else {
                if (last_exit) {
                    last_exit = nullptr;
                } else {
                    pos++;
                }
            }
        }

        cell.force_exit();

        if (cell.last_traversal().exited()) {
            // Record traversal for winding: get entry/exit y from the Traversal itself
            const auto& trav = cell.last_traversal();
            if (trav.traversed() && trav.coords().size() >= 2) {
                double t_entry_y = trav.coords().front().y;
                double t_exit_y  = trav.exit_coordinate().y;
                traversal_windings.push_back({row, col, t_entry_y, t_exit_y});
            }

            // Handle incomplete initial traversal (polygon start not on cell boundary)
            if (!cell.last_traversal().traversed()) {
                for (const auto& coord : cell.last_traversal().coords()) {
                    coords.push_back(coord);
                }
            }

            // Move to next cell based on exit side
            switch (cell.last_traversal().exit_side()) {
                case Side::TOP:    row--; break;
                case Side::BOTTOM: row++; break;
                case Side::LEFT:   col--; break;
                case Side::RIGHT:  col++; break;
                default: throw std::runtime_error("Invalid traversal exit side");
            }
        }
    }

    // ---- Extract coverage fractions and compute winding ----

    // For each cell that was visited, get coverage fraction
    for (auto& kv : cells) {
        size_t r = kv.first.first;
        size_t c = kv.first.second;

        float frac = static_cast<float>(kv.second->covered_fraction());
        if (frac == 0.0f) continue; // edge barely touched this cell

        // Map from infinite_extent grid coords to subgrid coords
        // In infinite_extent grid, the actual data starts at row=1, col=1
        // (padding=1 on each side)
        if (r < 1 || c < 1) continue;
        size_t sub_r = r - 1;
        size_t sub_c = c - 1;
        if (sub_r >= sub_rows || sub_c >= sub_cols) continue;

        // Find or create entry in row_data
        auto& row_vec = row_data[sub_r];
        // Look for existing entry for this column
        bool found = false;
        for (auto& rec : row_vec) {
            if (rec.col == static_cast<int>(col_off + sub_c)) {
                rec.coverage += coverage_factor * frac;
                found = true;
                break;
            }
        }
        if (!found) {
            row_vec.push_back({static_cast<int>(col_off + sub_c), coverage_factor * frac, 0});
        }
    }

    // Compute winding deltas from traversal records
    for (auto& tw : traversal_windings) {
        size_t r = tw.row;
        size_t c = tw.col;

        if (r < 1 || c < 1) continue;
        size_t sub_r = r - 1;
        size_t sub_c = c - 1;
        if (sub_r >= sub_rows || sub_c >= sub_cols) continue;

        // Row center y-coordinate
        // In the grid, row 0 (of the infinite_extent grid) is above ymax
        // Row 1 is the first actual row, its ymax = grid.ymax() - (but with padding offset)
        Box cell_box = grid_cell(grid, tw.row, tw.col);
        double y_mid = (cell_box.ymin + cell_box.ymax) / 2.0;

        // Does this traversal cross the row center?
        bool crosses = (tw.entry_y > y_mid && tw.exit_y < y_mid) ||
                       (tw.entry_y < y_mid && tw.exit_y > y_mid);

        if (!crosses) continue;

        // Direction: going downward (entry above mid, exit below) or upward?
        // For a CCW exterior ring:
        //   downward crossing → winding -1 (leaving interior on left)
        //   upward crossing   → winding +1 (entering interior on left)
        int delta;
        if (tw.entry_y > y_mid) {
            // Going downward
            delta = -1;
        } else {
            // Going upward
            delta = +1;
        }
        delta *= winding_factor;

        // Add winding delta to the corresponding row_data entry
        int full_col = static_cast<int>(col_off + sub_c);
        auto& row_vec = row_data[sub_r];
        bool found = false;
        for (auto& rec : row_vec) {
            if (rec.col == full_col) {
                rec.winding_delta += delta;
                found = true;
                break;
            }
        }
        if (!found) {
            // Coverage was zero but there's a winding crossing — this can happen
            // when an edge barely clips a cell corner
            row_vec.push_back({full_col, 0.0f, delta});
        }
    }
}

// Process a single GEOS geometry (polygon, multipolygon, etc)
static void process_geometry(
    GEOSContextHandle_t ctx,
    const GEOSGeometry* g,
    const exactextract::Grid<exactextract::bounded_extent>& full_grid,
    double dx, double dy,
    int poly_id,
    std::vector<GridRun>& all_runs,
    std::vector<GridEdge>& all_edges
) {
    using namespace exactextract;

    int type = GEOSGeomTypeId_r(ctx, g);

    if (type == GEOS_GEOMETRYCOLLECTION || type == GEOS_MULTIPOLYGON) {
        for (int i = 0; i < GEOSGetNumGeometries_r(ctx, g); i++) {
            process_geometry(ctx, GEOSGetGeometryN_r(ctx, g, i),
                             full_grid, dx, dy, poly_id, all_runs, all_edges);
        }
        return;
    }

    if (type != GEOS_POLYGON) {
        // For now, only handle polygons
        return;
    }

    // Get geometry bounding box
    auto component_boxes = geos_get_component_boxes(ctx, g);
    Box region = Box::make_empty();
    for (const auto& box : component_boxes) {
        if (!box.intersects(full_grid.extent())) continue;
        Box isect = full_grid.extent().intersection(box);
        if (region.empty()) {
            region = isect;
        } else if (!region.contains(isect)) {
            region = region.expand_to_include(isect);
        }
    }
    if (region.empty()) return;

    // Create subgrid with padding (infinite_extent)
    auto subgrid_bounded = full_grid.shrink_to_fit(region);
    auto subgrid = make_infinite(subgrid_bounded);

    if (subgrid.empty()) return;

    size_t sub_rows = subgrid.rows() - 2; // minus padding
    size_t sub_cols = subgrid.cols() - 2;

    // Compute offset of subgrid within full raster grid
    size_t row_off = static_cast<size_t>(
        std::round((full_grid.ymax() - subgrid_bounded.ymax()) / dy));
    size_t col_off = static_cast<size_t>(
        std::round((subgrid_bounded.xmin() - full_grid.xmin()) / dx));

    // Per-row boundary cell data (indexed by subgrid row)
    std::vector<std::vector<BoundaryCellRecord>> row_data(sub_rows);

    // Process exterior ring
    {
        const GEOSGeometry* ring = GEOSGetExteriorRing_r(ctx, g);
        const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq_r(ctx, ring);
        auto coords = read(ctx, seq);
        bool is_ccw = geos_is_ccw(ctx, seq);

        walk_ring(coords, is_ccw, true, subgrid, dy,
                  row_data, sub_rows, sub_cols, row_off, col_off);
    }

    // Process holes
    int n_holes = GEOSGetNumInteriorRings_r(ctx, g);
    for (int h = 0; h < n_holes; h++) {
        const GEOSGeometry* ring = GEOSGetInteriorRingN_r(ctx, g, h);
        const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq_r(ctx, ring);
        auto coords = read(ctx, seq);
        bool is_ccw = geos_is_ccw(ctx, seq);

        walk_ring(coords, is_ccw, false, subgrid, dy,
                  row_data, sub_rows, sub_cols, row_off, col_off);
    }

    // ---- Winding sweep: emit runs and edges per row ----

    float tol = 1e-6f;

    for (size_t sr = 0; sr < sub_rows; sr++) {
        auto& row_vec = row_data[sr];
        if (row_vec.empty()) continue;

        // Sort by column
        std::sort(row_vec.begin(), row_vec.end(),
            [](const BoundaryCellRecord& a, const BoundaryCellRecord& b) {
                return a.col < b.col;
            });

        // Merge entries for the same column (can happen with multi-ring cells)
        std::vector<BoundaryCellRecord> merged;
        for (auto& rec : row_vec) {
            if (!merged.empty() && merged.back().col == rec.col) {
                merged.back().coverage += rec.coverage;
                merged.back().winding_delta += rec.winding_delta;
            } else {
                merged.push_back(rec);
            }
        }

        // Sweep left to right with winding count
        int winding = 0;
        int prev_col = -2; // sentinel

        int full_row = static_cast<int>(row_off + sr) + 1; // 1-based

        for (auto& mc : merged) {
            // Emit interior run between previous boundary cell and this one
            if (winding != 0 && prev_col >= 0 && mc.col > prev_col + 1) {
                all_runs.push_back({
                    full_row,
                    prev_col + 1 + 1,   // 1-based, column after previous boundary
                    mc.col - 1 + 1,     // 1-based, column before this boundary
                    poly_id
                });
            }

            // Emit boundary cell
            float w = mc.coverage;
            if (w > tol && w < (1.0f - tol)) {
                // Partial coverage — emit as edge
                all_edges.push_back({full_row, mc.col + 1, w, poly_id});
            } else if (w >= (1.0f - tol)) {
                // Fully covered boundary cell — emit as single-cell run
                all_runs.push_back({full_row, mc.col + 1, mc.col + 1, poly_id});
            }
            // else: w <= tol means effectively outside, don't emit
            // (but still update winding below, since the edge still crosses)

            // Update winding count AFTER processing this cell
            winding += mc.winding_delta;
            prev_col = mc.col;
        }

        // Winding should return to 0 after all crossings on a row for a valid polygon
        // (don't enforce, just note if it's wrong)
    }
}


// ---- cpp11 entry point ----

[[cpp11::register]]
cpp11::writable::list cpp_scanline_burn(
    cpp11::list wkb_list,
    double xmin, double ymin, double xmax, double ymax,
    int ncol, int nrow
) {
    if (ncol <= 0 || nrow <= 0) {
        cpp11::stop("ncol and nrow must be positive");
    }
    if (xmax <= xmin || ymax <= ymin) {
        cpp11::stop("Invalid extent: xmax must be > xmin, ymax must be > ymin");
    }

    double dx = (xmax - xmin) / ncol;
    double dy = (ymax - ymin) / nrow;

    exactextract::Grid<exactextract::bounded_extent> full_grid(
        exactextract::Box(xmin, ymin, xmax, ymax), dx, dy);

    SLGEOSContextGuard geos_guard;
    GEOSContextHandle_t ctx = geos_guard.get();

    GEOSWKBReader* wkb_reader = GEOSWKBReader_create_r(ctx);
    if (!wkb_reader) cpp11::stop("Failed to create WKB reader");

    int n_geoms = wkb_list.size();

    std::vector<GridRun> all_runs;
    std::vector<GridEdge> all_edges;

    for (int k = 0; k < n_geoms; k++) {
        cpp11::raws wkb_raw(wkb_list[k]);
        int wkb_size = wkb_raw.size();
        if (wkb_size == 0) continue;

        const unsigned char* wkb_data = reinterpret_cast<const unsigned char*>(RAW(wkb_raw));

        SLGEOSGeomGuard geom(ctx,
            GEOSWKBReader_read_r(ctx, wkb_reader, wkb_data, static_cast<size_t>(wkb_size)));

        if (!geom.valid()) {
            cpp11::warning("Failed to parse WKB for geometry %d, skipping", k + 1);
            continue;
        }

        if (GEOSisEmpty_r(ctx, geom.get())) continue;

        try {
            int poly_id = k + 1; // 1-based
            process_geometry(ctx, geom.get(), full_grid, dx, dy,
                             poly_id, all_runs, all_edges);
        } catch (const std::exception& e) {
            cpp11::warning("Error processing geometry %d: %s, skipping", k + 1, e.what());
            continue;
        }
    }

    GEOSWKBReader_destroy_r(ctx, wkb_reader);

    // Build R data.frames (same format as cpp_burn_sparse)
    size_t n_runs = all_runs.size();
    cpp11::writable::integers runs_row(n_runs);
    cpp11::writable::integers runs_col_start(n_runs);
    cpp11::writable::integers runs_col_end(n_runs);
    cpp11::writable::integers runs_id(n_runs);

    for (size_t i = 0; i < n_runs; i++) {
        runs_row[i] = all_runs[i].row;
        runs_col_start[i] = all_runs[i].col_start;
        runs_col_end[i] = all_runs[i].col_end;
        runs_id[i] = all_runs[i].id;
    }

    cpp11::writable::list runs_df(4);
    runs_df[0] = static_cast<SEXP>(runs_row);
    runs_df[1] = static_cast<SEXP>(runs_col_start);
    runs_df[2] = static_cast<SEXP>(runs_col_end);
    runs_df[3] = static_cast<SEXP>(runs_id);
    runs_df.attr("names") = cpp11::writable::strings({"row", "col_start", "col_end", "id"});
    runs_df.attr("class") = "data.frame";
    runs_df.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n_runs)});

    size_t n_edges = all_edges.size();
    cpp11::writable::integers edges_row(n_edges);
    cpp11::writable::integers edges_col(n_edges);
    cpp11::writable::doubles edges_weight(n_edges);
    cpp11::writable::integers edges_id(n_edges);

    for (size_t i = 0; i < n_edges; i++) {
        edges_row[i] = all_edges[i].row;
        edges_col[i] = all_edges[i].col;
        edges_weight[i] = static_cast<double>(all_edges[i].weight);
        edges_id[i] = all_edges[i].id;
    }

    cpp11::writable::list edges_df(4);
    edges_df[0] = static_cast<SEXP>(edges_row);
    edges_df[1] = static_cast<SEXP>(edges_col);
    edges_df[2] = static_cast<SEXP>(edges_weight);
    edges_df[3] = static_cast<SEXP>(edges_id);
    edges_df.attr("names") = cpp11::writable::strings({"row", "col", "weight", "id"});
    edges_df.attr("class") = "data.frame";
    edges_df.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n_edges)});

    cpp11::writable::list result(2);
    result[0] = static_cast<SEXP>(runs_df);
    result[1] = static_cast<SEXP>(edges_df);
    result.attr("names") = cpp11::writable::strings({"runs", "edges"});

    return result;
}
