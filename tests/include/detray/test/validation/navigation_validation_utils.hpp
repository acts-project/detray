/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/barcode.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/tracks/free_track_parameters.hpp"
#include "detray/tracks/helix.hpp"

// Detray IO include(s)
#include "detray/io/csv/intersection2D.hpp"
#include "detray/io/csv/track_parameters.hpp"
#include "detray/io/utils/file_handle.hpp"

// Detray test include(s)
#include "detray/test/utils/inspectors.hpp"
#include "detray/test/validation/detector_scan_utils.hpp"
#include "detray/test/validation/material_validation_utils.hpp"
#include "detray/test/validation/step_tracer.hpp"

// Detray plugin include(s)
#include "detray/plugins/svgtools/styling/styling.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <iomanip>
#include <iterator>
#include <memory>
#include <ranges>
#include <sstream>

namespace detray::navigation_validator {

/// A functor to get the minimum distance to any surface boundary.
struct min_dist_to_boundary {

    template <typename mask_group_t, typename index_t, typename point_t>
    DETRAY_HOST_DEVICE inline auto operator()(const mask_group_t &mask_group,
                                              const index_t &index,
                                              const point_t &loc_p) const {

        return mask_group[index].min_dist_to_boundary(loc_p);
    }
};

/// B-field placeholder for straight-line navigation
struct empty_bfield {};

/// Struct that keeps track of the number of encountered/skipped surfaces
struct surface_stats {
    std::size_t n_portals{0u};
    std::size_t n_sensitives{0u};
    std::size_t n_passives{0u};

    /// The total number of skipped surfaces
    std::size_t n_total() const {
        return n_portals + n_sensitives + n_passives;
    }

    /// Count a surface depending on its type
    /// @returns true of the surface type was recognized
    /// @{
    template <typename sf_descriptor_t>
    bool count(const sf_descriptor_t &sf_desc) {
        return count(sf_desc.barcode());
    }

    bool count(geometry::barcode bcd) {
        switch (bcd.id()) {
            using enum surface_id;
            case e_portal: {
                n_portals++;
                return true;
            }
            case e_sensitive: {
                n_sensitives++;
                return true;
            }
            case e_passive: {
                n_passives++;
                return true;
            }
            default: {
                return false;
            }
        }
    }
    /// @}

    /// Count up
    surface_stats &operator+=(const surface_stats &other) {
        n_portals += other.n_portals;
        n_sensitives += other.n_sensitives;
        n_passives += other.n_passives;

        return *this;
    }
};

/// Statistics on tracks with holes and extra surfaces
struct track_statistics {
    std::size_t n_tracks{0u};
    std::size_t n_tracks_w_holes{0u};  //< Number of tracks with missing surface
    std::size_t n_tracks_w_extra{0u};  //< Number of tracks with extra surfaces
    std::size_t n_good_tracks{0u};     //< Number of tracks that match exactly
    std::size_t n_max_missed_per_trk{0u};  //< Max hole count
    std::size_t n_max_extra_per_trk{0u};   //< Max count of additional sf
};

/// Run the propagation and record test data along the way
template <typename stepper_t, typename... actor_ts, typename detector_t,
          typename bfield_t = empty_bfield>
inline auto record_propagation(
    const typename detector_t::geometry_context ctx,
    vecmem::memory_resource *host_mr, const detector_t &det,
    const propagation::config &cfg,
    const free_track_parameters<typename detector_t::algebra_type> &track,
    const bfield_t &bfield = {},
    const navigation::direction nav_dir = navigation::direction::e_forward) {

    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    /// Type that holds the intersection information
    using intersection_t =
        intersection2D<typename detector_t::surface_type, algebra_t, true>;

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
    using navigator_t = navigator<detector_t, navigation::default_cache_size,
                                  inspector_t, intersection_t>;

    // Propagator with pathlimit aborter and validation actors
    using step_tracer_t = step_tracer<algebra_t, dvector>;
    using pathlimit_aborter_t = pathlimit_aborter<scalar_t>;
    using material_tracer_t =
        material_validator::material_tracer<scalar_t, dvector>;
    using actor_chain_t = actor_chain<pathlimit_aborter_t, step_tracer_t,
                                      material_tracer_t, actor_ts...>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Propagator
    propagator_t prop{cfg};

    // Build actor states to collect data
    typename pathlimit_aborter_t::state pathlimit_aborter_state{
        cfg.stepping.path_limit};
    typename step_tracer_t::state step_tracer_state{*host_mr};
    typename material_tracer_t::state mat_tracer_state{*host_mr};

    // Extra actors
    dtuple<typename actor_ts::state...> state_tuple{};

    // Combine all actor states
    auto setup_actor_states =
        []<std::size_t... indices>(typename pathlimit_aborter_t::state & s1,
                                   typename step_tracer_t::state & s2,
                                   typename material_tracer_t::state & s3,
                                   dtuple<typename actor_ts::state...> & t,
                                   std::index_sequence<indices...> /*ids*/) {
        return detray::tie(s1, s2, s3, detail::get<indices>(t)...);
    };
    auto actor_states = setup_actor_states(
        pathlimit_aborter_state, step_tracer_state, mat_tracer_state,
        state_tuple, std::make_index_sequence<sizeof...(actor_ts)>{});

    std::unique_ptr<typename propagator_t::state> propagation{nullptr};
    if constexpr (std::is_same_v<bfield_t, empty_bfield>) {
        propagation =
            std::make_unique<typename propagator_t::state>(track, det, ctx);
    } else {
        propagation = std::make_unique<typename propagator_t::state>(
            track, bfield, det, ctx);
    }

    // Access to navigation information
    auto &nav_inspector = propagation->_navigation.inspector();
    auto &obj_tracer = nav_inspector.template get<object_tracer_t>();
    auto &nav_printer =
        nav_inspector.template get<navigation::print_inspector>();

    // Acces to the stepper information
    auto &step_printer = propagation->_stepping.inspector();

    // Find the end point and direction of the track (approximately)
    if (nav_dir == navigation::direction::e_backward) {
        using fw_propagator_t =
            propagator<stepper_t, navigator_t, actor_chain<>>;
        std::unique_ptr<typename fw_propagator_t::state> fw_propagation{
            nullptr};
        if constexpr (std::is_same_v<bfield_t, empty_bfield>) {
            fw_propagation = std::make_unique<typename fw_propagator_t::state>(
                track, det, ctx);
        } else {
            fw_propagation = std::make_unique<typename fw_propagator_t::state>(
                track, bfield, det, ctx);
        }
        fw_propagator_t{cfg}.propagate(*fw_propagation);
        propagation->_stepping() = fw_propagation->_stepping();
        propagation->_navigation.set_volume(
            fw_propagation->_navigation.volume());
    }

    // Run the propagation
    propagation->_navigation.set_direction(nav_dir);
    bool success = prop.propagate(*propagation, actor_states);

    return std::make_tuple(
        success, std::move(obj_tracer),
        std::move(step_tracer_state).release_step_data(),
        std::move(mat_tracer_state).release_material_record(),
        std::move(mat_tracer_state).release_material_steps(),
        std::move(nav_printer), std::move(step_printer));
}

/// Compare the recorded intersection trace in @param recorded_trace
/// to the truth trace in @param truth_trace
///
/// This function gathers and displays some statistics about missed and
/// additional surfaces in the recorded trace and provides detailed outputs
/// in case the traces do not match. If a hole is discovered in either trace,
/// a dummy record is added, so that both traces have the same length in the
/// end.
///
/// @param traj initial track parameters at the bginning of the track
/// @param trk_no number of the track in the sample
/// @param total_n_trks total number of tracks in the sample
/// @param debug_file where to output debugging information to
///
/// @returns the counts on missed and additional surfaces, matching errorsm
/// as well as the collections of the intersections that did not match
template <typename truth_trace_t, typename recorded_trace_t, typename traj_t>
auto compare_traces(truth_trace_t &truth_trace,
                    recorded_trace_t &recorded_trace, const traj_t &traj,
                    std::size_t trk_no, std::size_t total_n_trks,
                    std::fstream *debug_file = nullptr,
                    const bool fail_on_diff = true, const bool verbose = true) {

    using nav_record_t = typename recorded_trace_t::value_type;
    using truth_record_t = typename truth_trace_t::value_type;
    using intersection_t = typename truth_record_t::intersection_type;

    // Current number of entries to compare (may become more if both traces
    // missed several surfaces)
    std::size_t max_entries{
        math::max(recorded_trace.size(), truth_trace.size())};

    // Catch some debug output
    std::stringstream debug_stream;
    std::stringstream matching_stream;

    // Collect some statistics and additional data
    surface_stats missed_stats_nav{};
    surface_stats missed_stats_tr{};
    bool matching_traces{true};
    std::size_t n_errors{0u};
    std::vector<intersection_t> missed_intersections{};

    // An error might have occured occured
    auto handle_counting_error = [&matching_stream,
                                  &n_errors](const bool is_OK) {
        if (!is_OK) {
            matching_stream << "FATAL: Encountered surface of "
                               "unknown type in intersection "
                               "trace. Validate the geometry";
            ++n_errors;
        }
    };

    // Iterate until 'max_entries', because dummy records will be added to the
    // shorter trace
    for (long i = 0; i < static_cast<long>(max_entries); ++i) {

        // Check the records at the current index
        const auto idx{static_cast<std::size_t>(i)};

        // If only sensitives are collected and the recorded trace has one
        // entry less than the truth trace, the navigator missed the last
        // sensitive surface
        const bool nav_has_next = (idx < recorded_trace.size());
        detray::geometry::barcode nav_inters{};
        if (nav_has_next) {
            nav_inters = recorded_trace[idx].intersection.sf_desc.barcode();
        }

        const bool truth_has_next = (idx < truth_trace.size());
        detray::geometry::barcode truth_inters{};
        if (truth_has_next) {
            truth_inters = truth_trace[idx].intersection.sf_desc.barcode();
        }

        // Check if size of traces is still in sync and records match
        bool found_same_surfaces =
            (nav_has_next && truth_has_next && (nav_inters == truth_inters));

        matching_traces &= found_same_surfaces;

        if (!found_same_surfaces) {

            // Count the number of missed surfaces for this mismatch
            // Missed by navigator
            long missed_pt_nav{0};
            long missed_sn_nav{0};
            long missed_ps_nav{0};

            // Found in addition by navigation (missed by truth)
            long missed_pt_tr{0};
            long missed_sn_tr{0};
            long missed_ps_tr{0};

            // Intersection records at portal boundary might be flipped
            // (the portals overlap completely)
            auto is_swapped_portals = [&recorded_trace,
                                       &truth_trace](const long j) {
                const auto idx_j{static_cast<std::size_t>(j)};
                const std::size_t next_idx{idx_j + 1u};

                if (next_idx < truth_trace.size() &&
                    next_idx < recorded_trace.size()) {

                    const auto &current_nav_inters =
                        recorded_trace[idx_j].intersection.sf_desc.barcode();
                    const auto &current_truth_inters =
                        truth_trace[idx_j].intersection.sf_desc.barcode();

                    const auto &next_nav_inters =
                        recorded_trace[next_idx].intersection.sf_desc.barcode();
                    const auto &next_truth_inters =
                        truth_trace[next_idx].intersection.sf_desc.barcode();

                    return ((current_nav_inters == next_truth_inters) &&
                            (next_nav_inters == current_truth_inters));
                } else {
                    return false;
                }
            };

            // Count a missed surface by surface type
            auto count_one_missed = []<typename insers_t>(const insers_t &intr,
                                                          long int &missed_pt,
                                                          long int &missed_sn,
                                                          long int &missed_ps) {
                switch (intr.id()) {
                    using enum surface_id;
                    case e_portal: {
                        missed_pt++;
                        break;
                    }
                    case e_sensitive: {
                        missed_sn++;
                        break;
                    }
                    case e_passive: {
                        missed_ps++;
                        break;
                    }
                    default: {
                        throw std::runtime_error(
                            "Unkown surface type during counting");
                    }
                }
            };

            // Compare two traces and insert dummy records for any skipped cand.
            auto compare_and_equalize =
                [&i, &handle_counting_error, &is_swapped_portals,
                 &count_one_missed,
                 &missed_intersections]<typename trace_t,
                                        typename other_trace_t>(
                    trace_t &trace, typename trace_t::iterator last_missed_itr,
                    other_trace_t &other_trace, surface_stats &missed_stats,
                    long int &missed_pt, long int &missed_sn,
                    long int &missed_ps) {
                    // The navigator missed a(multiple) surface(s)
                    auto first_missed = std::begin(trace) + i;
                    const auto n_check{
                        std::distance(first_missed, last_missed_itr)};
                    assert(n_check > 0);

                    // Check and record surfaces that were missed: Insert dummy
                    // records until traces align again
                    for (long j = i; j < i + n_check; ++j) {
                        const auto &sfi =
                            trace[static_cast<std::size_t>(j)].intersection;

                        // Portals may be swapped and wrongfully included in the
                        // range of missed surfaces - skip them
                        if (sfi.sf_desc.is_portal() && is_swapped_portals(j)) {
                            ++j;
                            continue;
                        }
                        // Record the missed intersections fro later analysis
                        missed_intersections.push_back(sfi);
                        // Insert dummy record to match the truth trace size
                        using record_t = typename other_trace_t::value_type;
                        other_trace.insert(other_trace.begin() + i, record_t{});

                        // Count this missed intersection depending on sf. type
                        const bool valid{missed_stats.count(sfi.sf_desc)};
                        handle_counting_error(valid);

                        // Missed surfaces this time
                        count_one_missed(sfi.sf_desc, missed_pt, missed_sn,
                                         missed_ps);
                    }

                    assert(missed_pt >= 0);
                    assert(missed_sn >= 0);
                    assert(missed_ps >= 0);

                    const long n_missed{missed_pt + missed_sn + missed_ps};
                    // We landed here because something was missed
                    assert(n_missed > 0);

                    // Continue checking where trace might match again
                    i += (n_missed - 1);
                };

            // Match the barcodes to find how many surfaces were skipped
            //
            // If the current nav_inters can be found on the truth trace at a
            // later place, the navigator potentially missed the surfaces that
            // lie in between (except for swapped portals)
            auto search_nav_on_truth = [nav_inters](const truth_record_t &tr) {
                return tr.intersection.sf_desc.barcode() == nav_inters;
            };
            // As above, but this time check if the navigator found additional
            // surfaces
            auto search_truth_on_nav = [truth_inters](const nav_record_t &nr) {
                return nr.intersection.sf_desc.barcode() == truth_inters;
            };

            // Check if the portal order is swapped or the surface appears
            // later in the truth/navigation trace (this means one or
            // multiple surfaces were skipped respectively)
            if (is_swapped_portals(i)) {
                // Was not wrong after all
                matching_traces = true;
                // Have already checked the next record
                ++i;
            } else if (auto last_missed_by_nav = std::ranges::find_if(
                           std::ranges::begin(truth_trace) + i,
                           std::ranges::end(truth_trace), search_nav_on_truth);
                       last_missed_by_nav != std::end(truth_trace)) {

                // The navigator missed a(multiple) surface(s)
                compare_and_equalize(truth_trace, last_missed_by_nav,
                                     recorded_trace, missed_stats_nav,
                                     missed_pt_nav, missed_sn_nav,
                                     missed_ps_nav);

            } else if (auto last_missed_by_tr = std::ranges::find_if(
                           std::ranges::begin(recorded_trace) + i,
                           std::ranges::end(recorded_trace),
                           search_truth_on_nav);
                       last_missed_by_tr != std::end(recorded_trace)) {

                // The navigator found a(multiple) extra surface(s)
                compare_and_equalize(recorded_trace, last_missed_by_tr,
                                     truth_trace, missed_stats_tr, missed_pt_tr,
                                     missed_sn_tr, missed_ps_tr);

            } else if (!truth_has_next) {
                // The nav_inters could not be found on the truth trace, because
                // the truth trace does not have anymore records left to check:
                // The surface was missed on the truth side
                truth_trace.push_back(truth_record_t{});

                const bool valid{missed_stats_tr.count(nav_inters)};
                handle_counting_error(valid);
                if (valid) {
                    // Count to output error messages correctly
                    count_one_missed(nav_inters, missed_pt_tr, missed_sn_tr,
                                     missed_ps_tr);
                }
            } else if (!nav_has_next) {
                // The truth_inters could not be found on the recorded trace,
                // because the recorded trace does not have any records left
                // to check: The surface was missed by the navigator
                recorded_trace.push_back(nav_record_t{});

                const bool valid{missed_stats_nav.count(truth_inters)};
                handle_counting_error(valid);

                if (valid) {
                    // Count to output error messages correctly
                    count_one_missed(truth_inters, missed_pt_nav, missed_sn_nav,
                                     missed_ps_nav);
                }
            } else {
                // Both missed a surface at the same time, as neither record
                // can be found in each others traces
                bool valid{missed_stats_tr.count(nav_inters)};
                valid &= missed_stats_nav.count(truth_inters);
                handle_counting_error(valid);

                if (valid) {
                    // Count to output error messages correctly
                    count_one_missed(truth_inters, missed_pt_nav, missed_sn_nav,
                                     missed_ps_nav);

                    count_one_missed(nav_inters, missed_pt_tr, missed_sn_tr,
                                     missed_ps_tr);
                }
            }

            // Print error statements for the user
            auto print_err_extra_sf = [&matching_stream, max_entries, i](
                                          const std::string &sf_type,
                                          long n_sf) {
                matching_stream << "\nERROR: Detray navigator found " << n_sf
                                << " additional " << sf_type << "(s) at: " << i
                                << "/" << max_entries
                                << " (Inserted dummy record(s))";
            };

            auto print_err_missed = [&matching_stream, max_entries, i](
                                        const std::string &sf_type, long n_sf) {
                matching_stream << "\nERROR: Detray navigator missed " << n_sf
                                << " " << sf_type << "(s) at: " << i << "/"
                                << max_entries << ": "
                                << " (Inserted dummy record(s))";
            };

            if (missed_pt_tr > 0) {
                print_err_extra_sf("portal", missed_pt_tr);
            }
            if (missed_sn_tr > 0) {
                print_err_extra_sf("sensitive", missed_sn_tr);
            }
            if (missed_ps_tr > 0) {
                print_err_extra_sf("passive", missed_ps_tr);
            }

            if (missed_pt_nav > 0) {
                print_err_missed("portal", missed_pt_nav);
            }
            if (missed_sn_nav > 0) {
                print_err_missed("sensitive", missed_sn_nav);
            }
            if (missed_ps_nav > 0) {
                print_err_missed("passive", missed_ps_nav);
            }

            // Something must have been missed (unless it was just swapped
            // portals)
            assert(matching_traces ||
                   (missed_pt_tr + missed_sn_tr + missed_ps_tr + missed_pt_nav +
                        missed_sn_nav + missed_ps_nav >
                    0));
        }

        // Re-evaluate the size after dummy records were added
        max_entries = math::max(recorded_trace.size(), truth_trace.size());
    }

    matching_stream << "\n\nDetray navigator skipped "
                    << missed_stats_nav.n_total() << " surface(s) and found "
                    << missed_stats_tr.n_total() << " extra surface(s).";

    // Fill the debug stream with the final information from both traces
    for (std::size_t intr_idx = 0u; intr_idx < max_entries; ++intr_idx) {
        debug_stream << "-------Intersection ( " << intr_idx << " )\n";
        if (intr_idx < truth_trace.size()) {
            debug_stream << "\nReference: "
                         << truth_trace[intr_idx].intersection << ", vol id: "
                         << truth_trace[intr_idx].intersection.sf_desc.volume()
                         << std::endl;
        } else {
            debug_stream << "\nnReference: -" << std::endl;
        }
        if (intr_idx < recorded_trace.size()) {
            debug_stream << "\nDetray navigator:    "
                         << recorded_trace[intr_idx].intersection << std::endl
                         << std::endl;
        } else {
            debug_stream << "\nDetray navigator: -\n" << std::endl;
        }
    }

    const bool any_error{(missed_stats_nav.n_total() != 0u) ||
                         (missed_stats_tr.n_total() != 0u) || (n_errors != 0u)};

    // Fail the test with some extra information
    EXPECT_TRUE(!any_error || !fail_on_diff)
        << "\n--------\n"
        << "Track no. " << trk_no << "/" << total_n_trks << ":\n"
        << traj << matching_stream.str() << "\n--------";

    if (any_error && debug_file) {
        *debug_file << "\n>>>>>>>>>>>>>>>>>>\nFAILURE\n<<<<<<<<<<<<<<<<<<\n"
                    << "\nSUMMARY:\n--------\n"
                    << "Track no. " << trk_no << "/" << total_n_trks << ":\n"
                    << traj << matching_stream.str() << "\n--------\n"
                    << "\nFull Trace:\n\n"
                    << debug_stream.str();
    }

    // Multiple missed surfaces are a hint that something might be off with this
    // track
    if ((missed_stats_nav.n_total() > 1u) && verbose) {
        std::cout << "WARNING: Detray navigator skipped multiple surfaces: "
                  << missed_stats_nav.n_total() << "\n"
                  << std::endl;
    }
    if ((missed_stats_tr.n_total() > 1u) && verbose) {
        std::cout << "WARNING: Detray navigator found multiple extra surfaces: "
                  << missed_stats_tr.n_total() << "\n"
                  << std::endl;
    }

    // Unknown error occured during matching
    EXPECT_TRUE(n_errors == 0u)
        << "FATAL: Errors during matching: " << n_errors;

    // After inserting the placeholders, do a final check on the trace sizes
    const bool is_size{recorded_trace.size() == truth_trace.size()};
    EXPECT_TRUE(is_size)
        << "FATAL: Intersection traces have different number "
           "of surfaces after matching! Please check unmatched elements\n"
        << "Truth: " << truth_trace.size()
        << "\nNav. : " << recorded_trace.size() << "\n"
        << debug_stream.str();

    if (!is_size || (missed_stats_nav.n_total() != 0u) ||
        (missed_stats_tr.n_total() != 0u) || (n_errors != 0u)) {
        return std::make_tuple(false, missed_stats_nav, missed_stats_tr,
                               n_errors, missed_intersections);
    }

    // Make sure the failure was at least counted
    if (!matching_traces &&
        (missed_stats_nav.n_total() + missed_stats_tr.n_total() == 0)) {
        if (debug_file) {
            *debug_file << "\n>>>>>>>>>>>>>>>>>>\nFAILURE\n<<<<<<<<<<<<<<<<<<\n"
                        << "\nSUMMARY:\n--------\n"
                        << "Track no. " << trk_no << "/" << total_n_trks
                        << ":\n"
                        << traj << matching_stream.str() << "\n--------\n"
                        << "\nFull Trace:\n\n"
                        << debug_stream.str();
        }

        throw std::runtime_error(
            "Difference to truth trace was not counted correctly");
    }

    // Recount on the final traces, to be absolutely sure.
    if (!matching_traces) {
        std::size_t n_miss_truth{0};
        std::size_t n_miss_nav{0};
        for (std::size_t i = 0; i < truth_trace.size(); ++i) {
            const auto truth_desc{truth_trace[i].intersection.sf_desc};
            const auto nav_desc{recorded_trace[i].intersection.sf_desc};

            if (detail::is_invalid_value(truth_desc.volume())) {
                n_miss_truth++;
            } else if (detail::is_invalid_value(nav_desc.volume())) {
                n_miss_nav++;
            } else if (truth_desc.barcode() != nav_desc.barcode()) {
                n_miss_truth++;
                n_miss_nav++;
            }
        }

        if (n_miss_truth != missed_stats_tr.n_total()) {
            throw std::runtime_error(
                "Missed truth surfaces not counted correctly");
        }
        if (n_miss_nav != missed_stats_nav.n_total()) {
            throw std::runtime_error(
                "Missed navigation surfaces not counted correctly");
        }
    }

    const bool success{(!any_error || !fail_on_diff) && is_size};
    return std::make_tuple(success, missed_stats_nav, missed_stats_tr, n_errors,
                           missed_intersections);
}

/// Write the track positions of a trace @param intersection_traces to a csv
/// file to the path @param track_param_file_name
template <typename record_t>
auto write_tracks(const std::string &track_param_file_name,
                  const dvector<dvector<record_t>> &intersection_traces) {
    using algebra_t = typename record_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using track_param_t = free_track_parameters<algebra_t>;

    std::vector<std::vector<std::pair<scalar_t, track_param_t>>> track_params{};

    for (const auto &trace : intersection_traces) {

        track_params.push_back({});
        track_params.back().reserve(trace.size());

        for (const auto &record : trace) {
            track_params.back().emplace_back(
                record.charge,
                track_param_t{record.pos, 0.f, record.dir, record.charge});
        }
    }

    // Write to file
    io::csv::write_free_track_params(track_param_file_name, track_params);
}

/// Write the distance between the intersection and the surface boundaries in
/// @param missed_intersections to a csv file at the path @param file_name
template <typename detector_t, typename track_t, typename intersection_t>
auto write_dist_to_boundary(
    const detector_t &det, const typename detector_t::name_map &names,
    const std::string &file_name,
    const std::vector<std::pair<track_t, std::vector<intersection_t>>>
        &missed_intersections) {

    typename detector_t::geometry_context gctx{};

    // Write to csv file
    std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
    detray::io::file_handle dist_file{file_name, io_mode};
    *dist_file << "track_id,volume_id,volume_name,phi,eta,path,dist,inside_wo_"
                  "tol,sf_type"
               << std::endl;

    for (const auto &[i, entry] :
         detray::views::enumerate(missed_intersections)) {
        const auto &missed_inters_vec = entry.second;

        for (const auto &missed_sfi : missed_inters_vec) {

            const auto &track = entry.first;
            const auto sf =
                geometry::surface{det, missed_sfi.sf_desc.barcode()};
            const auto vol = tracking_volume{det, sf.volume()};

            const auto dist =
                sf.template visit_mask<min_dist_to_boundary>(missed_sfi.local);
            const auto glob_pos = sf.local_to_global(
                gctx, missed_sfi.local, track.dir(missed_sfi.path));

            *dist_file << i << "," << sf.volume() << ", " << vol.name(names)
                       << "," << vector::phi(glob_pos) << ", "
                       << vector::eta(glob_pos) << "," << missed_sfi.path
                       << ", " << dist << ", " << std::boolalpha
                       << sf.is_inside(missed_sfi.local, 0.f) << ", "
                       << static_cast<int>(sf.shape_id()) << std::endl;
        }
    }
}

/// Calculate and print the navigation efficiency
inline auto print_efficiency(std::size_t n_tracks,
                             const surface_stats &n_surfaces,
                             const surface_stats &n_miss_nav,
                             const surface_stats &n_miss_truth,
                             std::size_t n_fatal_error,
                             std::size_t n_matching_error) {
    // Column width in output
    constexpr int cw{20};

    // Print general information
    if (n_miss_nav.n_total() > 0u || n_miss_truth.n_total() > 0u ||
        n_fatal_error > 0u || n_matching_error > 0u) {

        std::cout
            << std::left << "-----------------------------------"
            << "Error Statistic:"
            << "\nTotal number of tracks: " << n_tracks
            << "\n\nTotal number of surfaces: " << n_surfaces.n_total()
            << std::setw(cw) << "\n      portals: " << n_surfaces.n_portals
            << std::setw(cw)
            << "\n      sensitives: " << n_surfaces.n_sensitives
            << std::setw(cw) << "\n      passives: " << n_surfaces.n_passives
            << "\n\n -> missed by navigator: " << n_miss_nav.n_total()
            << std::setw(cw) << "\n      portals: " << n_miss_nav.n_portals
            << std::setw(cw)
            << "\n      sensitives: " << n_miss_nav.n_sensitives
            << std::setw(cw) << "\n      passives: " << n_miss_nav.n_passives
            << "\n\n -> found in add. by navigator: " << n_miss_truth.n_total()
            << std::setw(cw) << "\n      portals: " << n_miss_truth.n_portals
            << std::setw(cw)
            << "\n      sensitives: " << n_miss_truth.n_sensitives
            << std::setw(cw) << "\n      passives: " << n_miss_truth.n_passives
            << "\n\nFatal propagation failures:   " << n_fatal_error
            << "\nErrors during truth matching: " << n_matching_error;
    } else {
        std::cout << "-----------------------------------\n"
                  << "Tested " << n_tracks << " tracks: OK\n"
                  << "total number of surfaces:         "
                  << n_surfaces.n_total();
    }

    assert(n_miss_nav.n_total() <= n_surfaces.n_total());
    assert(n_miss_nav.n_portals <= n_surfaces.n_portals);
    assert(n_miss_nav.n_sensitives <= n_surfaces.n_sensitives);
    assert(n_miss_nav.n_passives <= n_surfaces.n_passives);

    /// Print the surface finding efficiency per surface tpye
    auto print_eff = [&n_surfaces](const std::string &sf_type,
                                   const std::size_t n_sf,
                                   const std::size_t n_missed) {
        // How many significant digits to print
        const auto n_sig{2 + static_cast<int>(math::ceil(
                                 math::log10(n_surfaces.n_total())))};

        const auto k{static_cast<double>(n_sf - n_missed)};
        const auto n{static_cast<double>(n_sf)};

        // Estimate of the surface finding efficiency by the navigator
        const auto eff{k / n};

        // Variance
        // const double var_binomial{eff * (1. - eff) / n};
        const double var_bayesian{(k + 1.) * (k + 2.) / ((n + 2.) * (n + 3.)) -
                                  std::pow((k + 1.), 2) /
                                      std::pow((n + 2.), 2)};

        // In percent
        std::cout << "\n"
                  << sf_type << " finding eff.: " << std::fixed
                  << std::setprecision(n_sig) << 100. * eff << " \u00b1 "
                  << 100. * math::sqrt(var_bayesian) << "%";
    };

    std::cout << std::endl;
    if (n_surfaces.n_portals != 0u) {
        print_eff("Portal sf.", n_surfaces.n_portals, n_miss_nav.n_portals);
    }
    if (n_surfaces.n_sensitives != 0u) {
        print_eff("Sensitive sf.", n_surfaces.n_sensitives,
                  n_miss_nav.n_sensitives);
    }
    if (n_surfaces.n_passives != 0u) {
        print_eff("Passive sf.", n_surfaces.n_passives, n_miss_nav.n_passives);
    }

    std::cout << std::endl;
    if (n_surfaces.n_total() != 0u) {
        print_eff("Surface", n_surfaces.n_total(), n_miss_nav.n_total());
    } else {
        std::cout << "ERROR: No surfaces found in truth data!" << std::endl;
    }

    std::cout << "\n-----------------------------------\n" << std::endl;
}

/// Run the propagation and compare to an externally provided truth trace
///
template <typename stepper_t, typename... actor_ts, typename detector_t,
          typename field_view_t, typename intersection_t,
          concepts::algebra algebra_t>
auto compare_to_navigation(
    vecmem::host_memory_resource &host_mr, const detector_t &det,
    const typename detector_t::name_map &names,
    const typename detector_t::geometry_context ctx,
    const field_view_t field_view, const propagation::config &prop_cfg,
    std::vector<dvector<navigation::detail::candidate_record<intersection_t>>>
        &truth_traces,
    const std::vector<free_track_parameters<algebra_t>> &tracks,
    const navigation::direction nav_dir = navigation::direction::e_forward,
    const bool collect_sensitives_only = false, const bool fail_on_diff = true,
    const bool verbose = true) {

    // Style of svgs
    detray::svgtools::styling::style svg_style =
        detray::svgtools::styling::tableau_colorblind::style;

    // Write navigation and stepping debug info if a track fails
    std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
    detray::io::file_handle debug_file{"navigation_validation.txt", io_mode};

    // Collect some statistics
    std::size_t n_not_matching{0};
    std::size_t n_holes{0};
    std::size_t n_matching_error{0u};
    std::size_t n_fatal_error{0u};

    // Collect some statistics
    track_statistics trk_stats{};

    // Total number of encountered surfaces
    surface_stats n_surfaces{};
    // Missed by navigator
    surface_stats n_miss_nav{};
    // Missed by truth finder
    surface_stats n_miss_truth{};

    assert(truth_traces.size() == tracks.size());
    const std::size_t n_samples{truth_traces.size()};

    // Run the navigation on every truth particle
    for (std::size_t i = 0u; i < n_samples; ++i) {
        const auto &track = tracks[i];
        auto &truth_trace = truth_traces[i];

        assert(!truth_traces.empty());

        // Record the propagation through the geometry
        auto [success, obj_tracer, step_trace, mat_record, mat_trace,
              nav_printer, step_printer] =
            record_propagation<stepper_t, actor_ts...>(
                ctx, &host_mr, det, prop_cfg, track, field_view, nav_dir);

        // Fatal propagation error: Data unreliable
        if (!success) {
            std::cout << "ERROR: Propagation failure" << std::endl;

            *debug_file << "TEST TRACK " << i << ":\n\n"
                        << nav_printer.to_string() << step_printer.to_string();

            n_fatal_error++;
            n_not_matching++;
            continue;
        }

        // Compare to truth trace
        using obj_tracer_t = decltype(obj_tracer);
        using record_t = typename obj_tracer_t::candidate_record_t;

        detray::dvector<record_t> recorded_trace{};
        // Get only the sensitive surfaces
        if (collect_sensitives_only) {
            for (const auto &rec : obj_tracer.trace()) {
                if (rec.intersection.sf_desc.is_sensitive()) {
                    recorded_trace.push_back(rec);
                }
            }
        } else {
            std::ranges::copy(obj_tracer.trace(),
                              std::back_inserter(recorded_trace));
        }

        // Minimal number of incorrect traces (estimation)
        if (recorded_trace.size() != truth_trace.size()) {
            n_not_matching++;
        }
        if (recorded_trace.size() < truth_trace.size()) {
            n_holes++;
        }

        detray::detail::helix ideal_traj{track, field_view};
        auto [result, n_miss_trace_nav, n_miss_trace_truth, n_error_trace,
              missed_inters] =
            compare_traces(truth_trace, recorded_trace, ideal_traj, i,
                           n_samples, &(*debug_file), fail_on_diff, verbose);

        // Comparison failed
        if (!result) {
            // Write debug info to file
            *debug_file << "TEST TRACK " << i << ":\n\n"
                        << nav_printer.to_string() << step_printer.to_string();

            // In this test it is expected that the navigator finds
            // additional surfaces. Only dump svg if it missed one
            if (n_miss_trace_nav.n_total() != 0u || n_error_trace != 0u) {
                detray::detector_scanner::display_error(
                    ctx, det, names, "Navigation Check", ideal_traj,
                    truth_trace, svg_style, i, n_samples, recorded_trace);
            }
        }

        // After dummy records insertion, traces should have the same size
        EXPECT_EQ(truth_trace.size(), recorded_trace.size());

        // Add to global statistics
        n_miss_nav += n_miss_trace_nav;
        n_miss_truth += n_miss_trace_truth;
        n_matching_error += n_error_trace;

        // Count the number of different surface types on this trace
        for (std::size_t j = 0; j < truth_trace.size(); ++j) {
            n_surfaces.count(truth_trace[j].intersection.sf_desc);
        }

        // Did the navigation miss a sensitive surface? => count a hole
        if (n_miss_trace_nav.n_sensitives != 0u) {
            trk_stats.n_tracks_w_holes++;
        }
        if (n_error_trace == 0u && n_miss_trace_nav.n_total() == 0u &&
            n_miss_trace_truth.n_total() == 0u) {
            trk_stats.n_good_tracks++;
        }
        // Any additional surfaces found (might have material)
        if (n_miss_trace_truth.n_total() != 0u) {
            trk_stats.n_tracks_w_extra++;
        }

        trk_stats.n_max_missed_per_trk = math::max(
            trk_stats.n_max_missed_per_trk, n_miss_trace_nav.n_sensitives);
        trk_stats.n_max_extra_per_trk = math::max(trk_stats.n_max_extra_per_trk,
                                                  n_miss_trace_truth.n_total());

        if (n_error_trace != 0u) {
            throw std::runtime_error(
                "FATAL: Error during track comparison. Please check log "
                "files");
        }

        // Count the number of tracks that were tested
        trk_stats.n_tracks++;
    }

    const std::size_t n_tracks{trk_stats.n_tracks};
    const std::size_t n_tracks_w_holes{trk_stats.n_tracks_w_holes};
    const std::size_t n_tracks_w_extra{trk_stats.n_tracks_w_extra};
    const std::size_t n_good_tracks{trk_stats.n_good_tracks};
    const std::size_t n_max_holes_per_trk{trk_stats.n_max_missed_per_trk};
    const std::size_t n_max_extra_per_trk{trk_stats.n_max_extra_per_trk};

    // Self check
    if (n_not_matching > (n_tracks_w_holes + n_tracks_w_extra)) {
        throw std::runtime_error(
            "Number of tracks with mismatches underestimated!");
    }
    if (n_holes > n_tracks_w_holes) {
        throw std::runtime_error(
            "Number of tracks with holes is underestimated!");
    }

    // Calculate and display the surface finding efficiency
    print_efficiency(n_tracks, n_surfaces, n_miss_nav, n_miss_truth,
                     n_fatal_error, n_matching_error);

    // Column width in output
    constexpr int cw{35};

    std::cout << std::left << std::setw(cw)
              << "No. Tracks with holes: " << n_tracks_w_holes << "/"
              << n_tracks << " (" << 100. * n_tracks_w_holes / n_tracks << "%)"
              << std::endl;
    std::cout << std::left << std::setw(cw)
              << "No. Tracks with add. surfaces: " << n_tracks_w_extra << "/"
              << n_tracks << " (" << 100. * n_tracks_w_extra / n_tracks << "%)"
              << std::endl;
    std::cout << std::left << std::setw(cw)
              << "No. Good Tracks (exact match):  " << n_good_tracks << "/"
              << n_tracks << " (" << 100. * n_good_tracks / n_tracks << "%)\n"
              << std::endl;
    std::cout << std::left << std::setw(cw + 5)
              << "Max no. of holes per track: " << n_max_holes_per_trk
              << " (Mean: "
              << ((n_tracks_w_holes == 0)
                      ? "-"
                      : std::to_string(
                            static_cast<double>(n_miss_nav.n_sensitives) /
                            static_cast<double>(n_tracks_w_holes)))
              << ")" << std::endl;
    std::cout << std::left << std::setw(cw + 5)
              << "Max no. of add. surfaces per track: " << n_max_extra_per_trk
              << " (Mean: "
              << ((n_tracks_w_extra == 0)
                      ? "-"
                      : std::to_string(
                            static_cast<double>(n_miss_truth.n_total()) /
                            static_cast<double>(n_tracks_w_extra)))
              << ")" << std::endl;
    std::cout << "\n-----------------------------------" << std::endl;

    return std::tuple{trk_stats, n_surfaces, n_miss_nav, n_miss_truth};
}

}  // namespace detray::navigation_validator
