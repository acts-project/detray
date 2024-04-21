/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/tracks/free_track_parameters.hpp"
#include "detray/utils/inspectors.hpp"

// System include(s)
#include <algorithm>
#include <memory>
#include <sstream>

namespace detray::navigation_validator {

/// B-field placeholder for straight-line navigation
struct empty_bfield {};

/// Run the propagation and record test data along the way
template <typename stepper_t, typename detector_t,
          typename bfield_t = empty_bfield>
inline auto record_propagation(
    const typename detector_t::geometry_context, const detector_t &det,
    const propagation::config &cfg,
    const free_track_parameters<typename detector_t::algebra_type> &track,
    const bfield_t &bfield = {}) {

    using algebra_t = typename detector_t::algebra_type;

    /// Type that holds the intersection information
    using intersection_t =
        intersection2D<typename detector_t::surface_type, algebra_t>;

    /// Inspector that records all encountered surfaces
    using object_tracer_t =
        navigation::object_tracer<intersection_t, dvector,
                                  navigation::status::e_on_module,
                                  navigation::status::e_on_portal>;
    /// Inspector that prints the navigator state from within the
    /// navigator's method calls (cannot be done with an actor)
    using nav_print_inspector_t = navigation::print_inspector;
    /// Aggregation of multiple inspectors
    using inspector_t =
        aggregate_inspector<object_tracer_t, nav_print_inspector_t>;

    // Navigation with inspection
    using navigator_t = navigator<detector_t, inspector_t, intersection_t>;

    // Propagator with pathlimit aborter
    using actor_chain_t = actor_chain<dtuple, pathlimit_aborter>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Propagator
    propagator_t prop{cfg};

    // Build actor and propagator states
    pathlimit_aborter::state pathlimit_aborter_state{cfg.stepping.path_limit};
    auto actor_states = std::tie(pathlimit_aborter_state);

    std::unique_ptr<typename propagator_t::state> propagation{nullptr};
    if constexpr (std::is_same_v<bfield_t, empty_bfield>) {
        propagation =
            std::make_unique<typename propagator_t::state>(track, det);
    } else {
        propagation =
            std::make_unique<typename propagator_t::state>(track, bfield, det);
    }

    // Access to navigation information
    auto &nav_inspector = propagation->_navigation.inspector();
    auto &obj_tracer = nav_inspector.template get<object_tracer_t>();
    auto &nav_printer =
        nav_inspector.template get<navigation::print_inspector>();

    // Acces to the stepper information
    auto &step_printer = propagation->_stepping.inspector();

    // Run the propagation
    bool success = prop.propagate(*propagation, actor_states);

    return std::make_tuple(success, std::move(obj_tracer),
                           std::move(nav_printer), std::move(step_printer));
}

/// Compare the recorded intersection trace to the truth trace
template <typename truth_trace_t, typename recorded_trace_t, typename traj_t>
bool compare_traces(const truth_trace_t &truth_trace,
                    const recorded_trace_t &recorded_trace, const traj_t &traj,
                    std::size_t trk_no, std::size_t total_n_trks) {

    std::stringstream debug_stream;
    std::size_t n_inters_nav{recorded_trace.size()};
    std::size_t max_entries{math::max(n_inters_nav, truth_trace.size())};
    std::size_t min_entries{math::min(n_inters_nav, truth_trace.size())};

    // Fill the debug stream with the information from both traces
    for (std::size_t intr_idx = 0u; intr_idx < max_entries; ++intr_idx) {
        debug_stream << "-------Intersection ( " << intr_idx << " )\n";
        if (intr_idx < truth_trace.size()) {
            debug_stream << "\nparticle gun: "
                         << truth_trace[intr_idx].intersection << ", vol id: "
                         << truth_trace[intr_idx].intersection.sf_desc.volume()
                         << std::endl;
        } else {
            debug_stream << "\nparticle gun: -" << std::endl;
        }
        if (intr_idx < recorded_trace.size()) {
            debug_stream << "\nnavigator:    "
                         << recorded_trace[intr_idx].intersection << std::endl
                         << std::endl;
        } else {
            debug_stream << "\nnavigator: -\n" << std::endl;
        }
    }

    // Check every single recorded intersection
    for (std::size_t i = 0u; i < min_entries; ++i) {

        const auto &nav_inters =
            recorded_trace[i].intersection.sf_desc.barcode();
        const auto &ray_inters = truth_trace[i].intersection.sf_desc.barcode();

        const bool found_same_surfaces{nav_inters == ray_inters};

        if (not found_same_surfaces) {
            const auto &next_nav_inters =
                recorded_trace[i + 1u].intersection.sf_desc.barcode();
            const auto &next_ray_inters =
                truth_trace[i + 1u].intersection.sf_desc.barcode();

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
    const bool is_size{n_inters_nav == truth_trace.size()};
    EXPECT_TRUE(is_size) << "ERROR: Intersection traces found different number "
                            "of surfaces! Please check the last elements\n"
                         << debug_stream.str();
    if (not is_size) {
        return false;
    }

    return true;
}

/// Write the track positions of a trace @param intersection_traces to a csv
/// file to the path @param track_param_file_name
template <typename record_t>
auto write_tracks(const std::string &track_param_file_name,
                  const dvector<dvector<record_t>> &intersection_traces) {
    using track_param_t =
        free_track_parameters<typename record_t::algebra_type>;

    std::vector<std::vector<track_param_t>> track_params{};

    for (const auto &trace : intersection_traces) {

        track_params.push_back({});
        track_params.back().reserve(trace.size());

        for (const auto &record : trace) {
            track_params.back().push_back({record.pos, 0.f, record.dir, -1.f});
        }
    }

    // Write to file
    io::csv::write_free_track_params(track_param_file_name, track_params);
}

/// Calculate and print the navigation efficiency
/// @NOTE: WIP
inline auto print_efficiency(std::size_t n_tracks, std::size_t n_miss,
                             std::size_t n_fatal) {

    if (n_miss > 0u || n_fatal > 0u) {
        std::cout << "-----------------------------------"
                  << "Error Statistic:\n\n"
                  << "\n total: " << n_tracks << "\n (misses: " << n_miss
                  << ", fatal failures: " << n_fatal << ")\n"
                  << "-----------------------------------\n"
                  << std::endl;
    } else {
        std::cout << "-----------------------------------\n"
                  << "Tested " << n_tracks << " tracks: OK\n"
                  << "-----------------------------------\n"
                  << std::endl;
    }
}

}  // namespace detray::navigation_validator
