/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <algorithm>
#include <sstream>

namespace detray::detail {

template <typename inters_trace_t, typename object_tracer_t, typename traj_t>
bool compare_traces(const inters_trace_t &intersection_trace,
                    const object_tracer_t &obj_tracer, const traj_t &traj,
                    std::size_t trk_no, std::size_t total_n_trks) {

    std::stringstream debug_stream;
    std::size_t n_inters_nav{obj_tracer.object_trace.size()};
    std::size_t max_entries{math::max(n_inters_nav, intersection_trace.size())};
    std::size_t min_entries{math::min(n_inters_nav, intersection_trace.size())};

    // Fill the debug stream with the information from both traces
    for (std::size_t intr_idx = 0u; intr_idx < max_entries; ++intr_idx) {
        debug_stream << "-------Intersection ( " << intr_idx << " )\n";
        if (intr_idx < intersection_trace.size()) {
            debug_stream << "\nparticle gun: "
                         << intersection_trace[intr_idx].intersection
                         << ", vol id: " << intersection_trace[intr_idx].vol_idx
                         << std::endl;
        } else {
            debug_stream << "\nparticle gun: -" << std::endl;
        }
        if (intr_idx < obj_tracer.object_trace.size()) {
            debug_stream << "\nnavigator:    " << obj_tracer[intr_idx]
                         << std::endl
                         << std::endl;
        } else {
            debug_stream << "\nnavigator: -\n" << std::endl;
        }
    }

    // Check every single recorded intersection
    for (std::size_t i = 0u; i < min_entries; ++i) {

        const auto &nav_inters = obj_tracer[i].sf_desc.barcode();
        const auto &ray_inters =
            intersection_trace[i].intersection.sf_desc.barcode();

        const bool found_same_surfaces{nav_inters == ray_inters};

        if (not found_same_surfaces) {
            const auto &next_nav_inters = obj_tracer[i + 1u].sf_desc.barcode();
            const auto &next_ray_inters =
                intersection_trace[i + 1u].intersection.sf_desc.barcode();

            // Intersection record at portal bound might be flipped
            // (the portals overlap completely)
            if ((nav_inters == next_ray_inters) and
                (next_nav_inters == ray_inters)) {
                // Have already checked the next record
                ++i;
                continue;
            }
        }

        // Fail the test with some extra information
        EXPECT_TRUE(found_same_surfaces)
            << "\n>>>>>>>>>>>>>>>>>>\n"
            << "\nMismatch at intersection: " << i << "/" << n_inters_nav
            << " on track: " << trk_no << "/" << total_n_trks
            << "\n\nFailed for: " << traj << "\n<<<<<<<<<<<<<<<<<<\n"
            << "\nFull Trace:\n\n"
            << debug_stream.str();

        // Give the information about the failure to the caller
        if (not found_same_surfaces) {
            return false;
        }
    }

    // Do a final check on the trace sizes
    const bool is_size{n_inters_nav == intersection_trace.size()};
    EXPECT_TRUE(is_size) << "ERROR: Intersection traces found different number "
                            "of surfaces! Please check the last elements\n"
                         << debug_stream.str();
    if (not is_size) {
        return false;
    }

    return true;
}

}  // namespace detray::detail
