// analytical_coverage.h — analytical coverage fraction for single-traversal cells
//
// For a single edge traversal through a grid cell, the covered area (to the
// left of the traversal in CCW winding) is a simple polygon: the traversal
// path plus CCW cell boundary corners from exit to entry. This avoids the
// full left_hand_area chain-chasing algorithm in traversal_areas.cpp.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef DENSEBURN_ANALYTICAL_COVERAGE_H
#define DENSEBURN_ANALYTICAL_COVERAGE_H

#include <vector>
#include <cmath>
#include <cstddef>

#include "exactextract/box.h"
#include "exactextract/coordinate.h"
#include "exactextract/side.h"

namespace denseburn {

// Map Side to index in CCW perimeter order.
// CCW traversal of cell boundary: LEFT(up) → TOP(right) → RIGHT(down) → BOTTOM(left)
inline int side_ccw_index(exactextract::Side s) {
    switch (s) {
        case exactextract::Side::LEFT:   return 0;
        case exactextract::Side::TOP:    return 1;
        case exactextract::Side::RIGHT:  return 2;
        case exactextract::Side::BOTTOM: return 3;
        default: return -1;
    }
}

// Compute signed area of a closed polygon ring using the shoelace formula.
// Same formula as exactextract::area_signed (measures.cpp).
inline double polygon_signed_area(const std::vector<exactextract::Coordinate>& ring) {
    if (ring.size() < 3) return 0.0;

    double sum = 0.0;
    double x0 = ring[0].x;
    for (size_t i = 1; i < ring.size() - 1; i++) {
        double x = ring[i].x - x0;
        double y1 = ring[i + 1].y;
        double y2 = ring[i - 1].y;
        sum += x * (y2 - y1);
    }
    return sum / 2.0;
}

// Compute exact coverage fraction for a single traversal through a cell.
//
// The traversal enters from entry_side at coords.front() and exits through
// exit_side at coords.back(). Intermediate points (polygon vertices within
// the cell) may appear between entry and exit.
//
// For a CCW-oriented ring, the covered area is to the LEFT of the traversal
// direction. This equals the area of the polygon formed by:
//   1. The traversal path (entry → ... → exit)
//   2. CCW cell boundary from exit back to entry (inserting corner points)
//
// Returns coverage fraction in [0, 1].
inline double analytical_covered_fraction(
    const exactextract::Box& box,
    const std::vector<exactextract::Coordinate>& coords,
    exactextract::Side entry_side,
    exactextract::Side exit_side
) {
    using exactextract::Coordinate;

    // Build the left-hand polygon
    std::vector<Coordinate> polygon;
    polygon.reserve(coords.size() + 5); // traversal + up to 4 corners + close

    // 1. Traversal coordinates
    polygon.insert(polygon.end(), coords.begin(), coords.end());

    // 2. CCW boundary from exit to entry, inserting corners
    //
    // Corners in CCW order, indexed by the side they follow:
    //   LEFT(0)  → TL (xmin, ymax)
    //   TOP(1)   → TR (xmax, ymax)
    //   RIGHT(2) → BR (xmax, ymin)
    //   BOTTOM(3)→ BL (xmin, ymin)
    const Coordinate corners[4] = {
        {box.xmin, box.ymax},  // TL: CCW end of LEFT
        {box.xmax, box.ymax},  // TR: CCW end of TOP
        {box.xmax, box.ymin},  // BR: CCW end of RIGHT
        {box.xmin, box.ymin},  // BL: CCW end of BOTTOM
    };

    int exit_idx = side_ccw_index(exit_side);
    int entry_idx = side_ccw_index(entry_side);

    // Walk CCW from exit side to entry side, adding corner at each transition
    int i = exit_idx;
    while (i != entry_idx) {
        polygon.push_back(corners[i]);
        i = (i + 1) % 4;
    }

    // Close the polygon
    polygon.push_back(polygon.front());

    // Compute area
    double cell_area = box.area();
    if (cell_area <= 0.0) return 0.0;

    return std::abs(polygon_signed_area(polygon)) / cell_area;
}

// Compute coverage fraction for a closed ring entirely within one cell.
// This handles the case where a small polygon fits inside a single grid cell
// (no entry/exit sides — the ring is self-contained).
inline double closed_ring_covered_fraction(
    const exactextract::Box& box,
    const std::vector<exactextract::Coordinate>& ring_coords
) {
    double cell_area = box.area();
    if (cell_area <= 0.0) return 0.0;
    return std::abs(polygon_signed_area(ring_coords)) / cell_area;
}

} // namespace denseburn

#endif // DENSEBURN_ANALYTICAL_COVERAGE_H
