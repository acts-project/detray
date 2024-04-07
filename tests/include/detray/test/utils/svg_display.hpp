/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/utils/create_path.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// System include(s)
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_set>

namespace detray::detail {

/// Find the unique volume indices that the trajectory crossed
/// - intersection record
template <typename record_t>
std::unordered_set<dindex> get_volume_indices(
    const std::vector<record_t> &intersection_record) {

    std::unordered_set<dindex> volumes{};
    volumes.reserve(intersection_record.size());
    for (const auto &single_ir : intersection_record) {
        volumes.insert(single_ir.vol_idx);
    }

    return volumes;
}

/// Find the unique volume indices that the trajectory crossed
/// - intersection collection
template <typename surface_t, typename algebra_t>
std::unordered_set<dindex> get_volume_indices(
    const dvector<detray::intersection2D<surface_t, algebra_t>>
        &intersections) {

    std::unordered_set<dindex> volumes{};
    volumes.reserve(intersections.size());
    for (const auto &intr : intersections) {
        volumes.insert(intr.sf_desc.volume());
    }

    return volumes;
}

/// @returns the svg of the intersections (truth and track) and the trajectory
template <typename detector_t, typename record_t, class traj_t, typename view_t>
auto draw_intersection_and_traj_svg(
    const typename detector_t::geometry_context gctx,
    detray::svgtools::illustrator<detector_t> &il,
    const std::vector<record_t> &intersections_truth, const traj_t &traj,
    const std::string &traj_name,
    const dvector<typename record_t::intersection_type> &intersections,
    const view_t &view) {

    auto svg_traj = il.draw_intersections(
        "truth_intersections", intersections_truth, traj.dir(), view, gctx);

    if (not intersections.empty()) {
        svg_traj.add_object(il.draw_intersections_and_trajectory(
            traj_name, intersections, traj, view,
            intersections_truth.back().intersection.path, gctx));
    } else {
        svg_traj.add_object(il.draw_trajectory(
            traj_name, traj, intersections_truth.back().intersection.path,
            view));
    }

    return svg_traj;
}

/// Display the geometry, intersection and track data via @c svgtools
template <typename detector_t, typename record_t, class traj_t>
inline void svg_display(
    const typename detector_t::geometry_context gctx,
    detray::svgtools::illustrator<detector_t> &il,
    const std::vector<record_t> &intersections_truth, const traj_t &traj,
    const std::string &traj_name,
    const std::string &outfile = "detector_display",
    const dvector<typename record_t::intersection_type> &intersections = {},
    const std::string &outdir = "./plots/") {

    // Gather all volumes that need to be displayed
    auto volumes = get_volume_indices(intersections_truth);
    if (not intersections.empty()) {
        const auto more_volumes = get_volume_indices(intersections_truth);
        volumes.insert(more_volumes.begin(), more_volumes.end());
    }

    // General options
    auto path = detray::io::create_path(outdir);

    actsvg::style::stroke stroke_black = actsvg::style::stroke();

    // x-y axis.
    auto xy_axis = actsvg::draw::x_y_axes("axes", {-1100, 1100}, {-1100, 1100},
                                          stroke_black, "x", "y");
    // z-r axis.
    auto zr_axis = actsvg::draw::x_y_axes("axes", {-3100, 3100}, {-5, 1100},
                                          stroke_black, "z", "r");
    // Creating the views.
    const actsvg::views::x_y xy;
    const actsvg::views::z_r zr;

    // xy - view
    auto svg_traj = draw_intersection_and_traj_svg(
        gctx, il, intersections_truth, traj, traj_name, intersections, xy);

    const auto [vol_xy_svg, _] = il.draw_volumes(volumes, xy, gctx);
    detray::svgtools::write_svg(
        path / (outfile + "_" + vol_xy_svg._id + "_" + traj_name),
        {xy_axis, vol_xy_svg, svg_traj});

    // zr - view
    svg_traj = draw_intersection_and_traj_svg(
        gctx, il, intersections_truth, traj, traj_name, intersections, zr);

    const auto vol_zr_svg = il.draw_detector(zr, gctx);
    detray::svgtools::write_svg(
        path / (outfile + "_" + vol_zr_svg._id + "_" + traj_name),
        {zr_axis, vol_zr_svg, svg_traj});

    std::cout << "INFO: Wrote debugging data in: " << path << "\n" << std::endl;
}

}  // namespace detray::detail
