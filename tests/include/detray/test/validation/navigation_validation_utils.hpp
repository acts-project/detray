/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/track_parametrization.hpp"
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
#include "detray/test/common/navigation_validation_config.hpp"
#include "detray/test/utils/inspectors.hpp"
#include "detray/test/validation/detector_scan_utils.hpp"
#include "detray/test/validation/material_validation_utils.hpp"
#include "detray/test/validation/step_tracer.hpp"

// Detray plugin include(s)
#include "detray/plugins/svgtools/styling/styling.hpp"

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
    bool count(const sf_descriptor_t &sf_desc, const bool verbose = true) {
        return count(sf_desc.barcode(), verbose);
    }

    bool count(geometry::barcode bcd, const bool verbose = true) {
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
                if (verbose) {
                    std::cout << "WARNING: Surface type unknown " << bcd
                              << std::endl;
                }
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
    const pdg_particle<typename detector_t::scalar_type> ptc_hypo =
        muon<typename detector_t::scalar_type>(),
    const bfield_t &bfield = {},
    const navigation::direction nav_dir = navigation::direction::e_forward,
    typename actor_chain<actor_ts...>::state_ref_tuple state_tuple = {},
    const std::array<dscalar<typename detector_t::algebra_type>, e_bound_size>
        &stddevs = {}) {

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
    using actor_chain_t = actor_chain<pathlimit_aborter_t, actor_ts...,
                                      step_tracer_t, material_tracer_t>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Propagator
    propagator_t prop{cfg};

    // Build actor states to collect data
    typename pathlimit_aborter_t::state pathlimit_aborter_state{
        cfg.stepping.path_limit};
    typename step_tracer_t::state step_tracer_state{*host_mr};
    typename material_tracer_t::state mat_tracer_state{*host_mr};

    // Combine all actor states
    auto setup_actor_states = []<std::size_t... indices>(
        typename pathlimit_aborter_t::state & s1,
        typename step_tracer_t::state & s2,
        typename material_tracer_t::state & s3,
        typename actor_chain<actor_ts...>::state_ref_tuple & t,
        std::index_sequence<indices...> /*ids*/) {
        return detray::tie(s1, detail::get<indices>(t)..., s2, s3);
    };
    constexpr auto idx_seq{std::make_index_sequence<detail::tuple_size_v<
        typename actor_chain<actor_ts...>::state_ref_tuple>>{}};

    auto actor_states =
        setup_actor_states(pathlimit_aborter_state, step_tracer_state,
                           mat_tracer_state, state_tuple, idx_seq);

    std::unique_ptr<typename propagator_t::state> propagation{nullptr};
    if constexpr (std::is_same_v<bfield_t, empty_bfield>) {
        propagation =
            std::make_unique<typename propagator_t::state>(track, det, ctx);
    } else {
        propagation = std::make_unique<typename propagator_t::state>(
            track, bfield, det, ctx);
    }

    // Set the initial covariances
    if (!stddevs.empty()) {
        std::random_device rd{};
        std::mt19937 generator{rd()};

        auto &bound_param = propagation->_stepping.bound_params();

        for (std::size_t i = 0u; i < e_bound_size; i++) {

            if (stddevs[i] != scalar_t{0}) {
                bound_param[i] = std::normal_distribution<scalar_t>(
                    bound_param[i], stddevs[i])(generator);
            }

            getter::element(bound_param.covariance(), i, i) =
                stddevs[i] * stddevs[i];
        }
    }

    // Access to navigation information
    auto &nav_inspector = propagation->_navigation.inspector();
    auto &obj_tracer = nav_inspector.template get<object_tracer_t>();
    auto &nav_printer =
        nav_inspector.template get<navigation::print_inspector>();

    // Acces to the stepper information
    auto &step_printer = propagation->_stepping.inspector();

    // Find the end point and direction of the track (approximately) for
    // backward propagation
    if (nav_dir == navigation::direction::e_backward) {
        using fw_propagator_t =
            propagator<stepper_t, navigator_t,
                       actor_chain<pathlimit_aborter_t, actor_ts...>>;
        std::unique_ptr<typename fw_propagator_t::state> fw_propagation{
            nullptr};
        if constexpr (std::is_same_v<bfield_t, empty_bfield>) {
            fw_propagation = std::make_unique<typename fw_propagator_t::state>(
                track, det, ctx);
        } else {
            fw_propagation = std::make_unique<typename fw_propagator_t::state>(
                track, bfield, det, ctx);
        }

        // Make a deep copy of states for forward propagation, but omit the
        // expensive tracers for the forward pass
        auto copy_actor_states = []<std::size_t... indices>(
            typename pathlimit_aborter_t::state & s1,
            typename actor_chain<actor_ts...>::state_ref_tuple & t,
            std::index_sequence<indices...> /*ids*/) {
            using fw_state_tuple_t =
                typename actor_chain<pathlimit_aborter_t,
                                     actor_ts...>::state_tuple;
            return fw_state_tuple_t{s1, detail::get<indices>(t)...};
        };

        // Setup forward actor states
        auto setup_fw_actor_states = []<std::size_t... indices>(
            typename actor_chain<pathlimit_aborter_t,
                                 actor_ts...>::state_tuple &
                t,
            std::index_sequence<indices...> /*ids*/) {
            return detray::tie(detail::get<indices>(t)...);
        };

        auto fw_state_tuple =
            copy_actor_states(pathlimit_aborter_state, state_tuple, idx_seq);
        constexpr auto fw_idx_seq{std::make_index_sequence<
            detail::tuple_size_v<decltype(fw_state_tuple)>>{}};
        auto fw_actor_states =
            setup_fw_actor_states(fw_state_tuple, fw_idx_seq);

        // Perform forward propagation
        fw_propagation->set_particle(
            update_particle_hypothesis(ptc_hypo, track));
        fw_propagation->_stepping.bound_params() =
            propagation->_stepping.bound_params();

        const bool fw_success =
            fw_propagator_t{cfg}.propagate(*fw_propagation, fw_actor_states);

        if (!fw_success) {
            std::cout << "ERROR: Could not propagate to end of track to "
                         "prepare backward propagation"
                      << std::endl;
        }

        // USe the result to set up main propagation run
        propagation->_stepping() = fw_propagation->_stepping();
        propagation->_stepping.bound_params() =
            fw_propagation->_stepping.bound_params();
        propagation->_navigation.set_volume(
            fw_propagation->_navigation.volume());
    }

    // Run the propagation
    propagation->_navigation.set_direction(nav_dir);
    propagation->set_particle(update_particle_hypothesis(ptc_hypo, track));
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
/// @param debug_file where to output debugging information to
///
/// @returns the counts on missed and additional surfaces, matching errorsm
/// as well as the collections of the intersections that did not match
template <typename truth_trace_t, typename recorded_trace_t, typename traj_t>
auto compare_traces(const detray::test::navigation_validation_config &cfg,
                    truth_trace_t &truth_trace,
                    recorded_trace_t &recorded_trace, const traj_t &traj,
                    std::size_t trk_no, std::fstream *debug_file = nullptr) {

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

    // Missed surfaces or errors during matching
    const bool any_error{(missed_stats_nav.n_total() != 0u) ||
                         (missed_stats_tr.n_total() != 0u) || (n_errors != 0u)};

    // Print debugging information if anything went wrong
    if (any_error && cfg.verbose()) {
        std::cout << "\n--------\n"
                  << "Track no. " << trk_no << "/" << cfg.n_tracks() << ":\n"
                  << traj << matching_stream.str() << "\n--------";
        if (debug_file) {
            *debug_file << "\n>>>>>>>>>>>>>>>>>>\nFAILURE\n<<<<<<<<<<<<<<<<<<\n"
                        << "\nSUMMARY:\n--------\n"
                        << "Track no. " << trk_no << "/" << cfg.n_tracks()
                        << ":\n"
                        << traj << matching_stream.str() << "\n--------\n"
                        << "\nFull Trace:\n\n"
                        << debug_stream.str();
        }
    }

    // Unknown error occured during matching
    if (n_errors != 0u) {
        std::stringstream err_str{};
        err_str << "FATAL: Errors during matching: " << n_errors << "\n"
                << matching_stream.str();
        throw std::runtime_error(err_str.str());
    }

    // After inserting the placeholders, do a final check on the trace sizes
    const bool is_size{recorded_trace.size() == truth_trace.size()};
    if (!is_size) {
        std::stringstream err_str{};
        err_str
            << "FATAL: Intersection traces have different number "
               "of surfaces after matching! Please check unmatched elements\n"
            << "Truth: " << truth_trace.size()
            << "\nNav. : " << recorded_trace.size() << "\n"
            << debug_stream.str();
        throw std::runtime_error(err_str.str());
    }

    // Multiple missed surfaces are a hint that something might be off with this
    // track
    if ((missed_stats_nav.n_total() > 1u) && cfg.verbose()) {
        std::cout << "WARNING: Detray navigator skipped multiple surfaces: "
                  << missed_stats_nav.n_total() << "\n"
                  << std::endl;
    }
    if ((missed_stats_tr.n_total() > 1u) && cfg.verbose()) {
        std::cout << "WARNING: Detray navigator found multiple extra surfaces: "
                  << missed_stats_tr.n_total() << "\n"
                  << std::endl;
    }

    // Make sure the mismatches were correctly counted
    if (!matching_traces) {
        bool is_counting_error{false};
        // No mismatch counted for traces that don't match: Error
        if (missed_stats_nav.n_total() + missed_stats_tr.n_total() == 0) {
            is_counting_error = true;

            const std::string msg{
                "Difference to truth trace was not counted correctly"};
            if (cfg.fail_on_diff()) {
                throw std::runtime_error(msg);
            } else {
                std::cout << "ERROR: " << msg << std::endl;
            }
        }

        // Recount on the final traces, to be absolutely sure.
        std::size_t n_miss_truth{0u};
        std::size_t n_miss_nav{0u};
        for (std::size_t i = 0u; i < truth_trace.size(); ++i) {
            const auto truth_desc{truth_trace[i].intersection.sf_desc};
            const auto nav_desc{recorded_trace[i].intersection.sf_desc};

            if (truth_desc.barcode().is_invalid()) {
                n_miss_truth++;
            } else if (nav_desc.barcode().is_invalid()) {
                n_miss_nav++;
            } else if (truth_desc.barcode() != nav_desc.barcode()) {
                n_miss_truth++;
                n_miss_nav++;
            }
        }

        if (n_miss_truth != missed_stats_tr.n_total()) {
            is_counting_error = true;
            const std::string msg{
                "Missed truth surfaces not counted correctly"};
            if (cfg.fail_on_diff()) {
                throw std::runtime_error(msg);
            } else {
                std::cout << "ERROR: " << msg << std::endl;
            }
        }
        if (n_miss_nav != missed_stats_nav.n_total()) {
            is_counting_error = true;
            const std::string msg{
                "Missed navigation surfaces not counted correctly"};
            if (cfg.fail_on_diff()) {
                throw std::runtime_error(msg);
            } else {
                std::cout << "ERROR: " << msg << std::endl;
            }
        }

        if (is_counting_error) {
            n_errors++;

            if (debug_file) {
                *debug_file << "\n>>>>>>>>>>>>>>>>>>\nCOUNTING "
                               "ERROR\n<<<<<<<<<<<<<<<<<<\n"
                            << "\nSUMMARY:\n--------\n"
                            << "Track no. " << trk_no << "/" << cfg.n_tracks()
                            << ":\n"
                            << traj << matching_stream.str() << "\n--------\n"
                            << "\nFull Trace:\n\n"
                            << debug_stream.str();
            }
        }
    }

    const bool success{!any_error && is_size};
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

/// Run the propagation and compare to an externally provided truth trace.
///
/// The resulting statistics on tracks with missed surfaces etc. will be
/// printed to stdout. For any track with mismatches in the surface trace, a
/// detailed report will be dumped into a debug file and an svg showing the
/// track and truth trace in the detector geometry will be provided in ./plots
///
/// @tparam stepper_t the type of stepper to use (e.g. straight line vs. RKN)
/// @tparam actor_ts types of additional actors (e.g. parameter transport)
///
/// @param cfg navigation validation config
/// @param host_mr host memory resource
/// @param det the detector
/// @param names voluem names for the detector
/// @param ctx the geometry context
/// @param field_view magnetic field view
/// @param prop_cfg the propagation configuration
/// @param truth_traces coll. of truth traces (one per track)
/// @param tracks coll. of tracks
/// @param state_tuples coll. of actor state tuples for actor_ts (one per track)
/// @param stddevs_per_track coll. standard dev. for bound parameter smearing
///
/// @returns tuple of track statistics: stats of tracks (holes etc.), stats of
/// encountered surfaces, stats of missed surfaces for navigation, stats of
/// missed surfaces for truth traces, recorded step traces, recorded material
/// traces and intergrated material
template <typename stepper_t, typename... actor_ts, typename detector_t,
          typename field_view_t, typename intersection_t,
          concepts::algebra algebra_t>
auto compare_to_navigation(
    const detray::test::navigation_validation_config &cfg,
    vecmem::host_memory_resource &host_mr, const detector_t &det,
    const typename detector_t::name_map &names,
    const typename detector_t::geometry_context ctx,
    const field_view_t field_view, const propagation::config &prop_cfg,
    std::vector<dvector<navigation::detail::candidate_record<intersection_t>>>
        &truth_traces,
    const std::vector<free_track_parameters<algebra_t>> &tracks,
    dvector<typename actor_chain<actor_ts...>::state_ref_tuple> state_tuples =
        {{}},
    const dvector<std::array<dscalar<algebra_t>, e_bound_size>>
        &stddevs_per_track = {{0.f}}) {

    using scalar_t = dscalar<algebra_t>;

    assert(truth_traces.size() == tracks.size());
    assert(!state_tuples.empty());
    assert(!stddevs_per_track.empty());

    // Write navigation and stepping debug info if a track fails
    std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
    detray::io::file_handle debug_file{
        "navigation_validation_" + cfg.name() + ".txt", io_mode};

    const std::size_t n_samples{truth_traces.size()};

    // Collect some statistics for consistency checking
    std::size_t n_holes{0};
    std::size_t n_matching_error{0u};
    std::size_t n_fatal_error{0u};

    // Collect global track statistics
    track_statistics trk_stats{};

    // Total number of encountered surfaces
    surface_stats n_surfaces{};
    // Missed by navigator
    surface_stats n_miss_nav{};
    // Missed in truth trace
    surface_stats n_miss_truth{};

    // Collect step and material traces for all tracks
    std::vector<dvector<detail::step_data<algebra_t>>> step_traces{};
    std::vector<material_validator::material_record<scalar_t>> mat_records{};
    std::vector<dvector<material_validator::material_params<scalar_t>>>
        mat_traces{};

    mat_records.reserve(n_samples);
    mat_traces.reserve(n_samples);
    step_traces.reserve(n_samples);

    // Run the navigation on every truth trace
    for (std::size_t i = 0u; i < n_samples; ++i) {
        const auto &track = tracks.at(i);
        auto &truth_trace = truth_traces.at(i);
        auto state_tuple =
            state_tuples.size() > i ? state_tuples.at(i) : state_tuples.at(0);
        const auto &stddevs = stddevs_per_track.size() > i
                                  ? stddevs_per_track.at(i)
                                  : std::array<scalar_t, e_bound_size>{0.f};

        assert(!truth_traces.empty());

        // Record the propagation through the geometry
        auto [success, obj_tracer, step_trace, mat_record, mat_trace,
              nav_printer, step_printer] =
            record_propagation<stepper_t, actor_ts...>(
                ctx, &host_mr, det, prop_cfg, track, cfg.ptc_hypothesis(),
                field_view, cfg.navigation_direction(), state_tuple, stddevs);

        // Fatal propagation error: Data unreliable
        if (!success) {
            std::cout << "ERROR: Propagation failure" << std::endl;

            *debug_file << "ERROR: Propagation failure:\n"
                        << "TEST TRACK " << i;

            n_fatal_error++;
        }

        // Save material for later comparison
        step_traces.push_back(std::move(step_trace));
        mat_records.push_back(std::move(mat_record));
        mat_traces.push_back(mat_trace);

        // Compare to truth trace
        using obj_tracer_t = decltype(obj_tracer);
        using record_t = typename obj_tracer_t::candidate_record_t;

        detray::dvector<record_t> recorded_trace{};
        recorded_trace.reserve(obj_tracer.trace().size());
        // Get only the sensitive surfaces
        if (cfg.collect_sensitives_only()) {
            for (const auto &rec : obj_tracer.trace()) {
                if (rec.intersection.sf_desc.is_sensitive()) {
                    recorded_trace.push_back(rec);
                }
            }
        } else {
            std::ranges::copy(obj_tracer.trace(),
                              std::back_inserter(recorded_trace));
        }

        // The recorded trace missed somethin, make an estimate of min. number
        // of missed surfaces
        if (recorded_trace.size() < truth_trace.size()) {
            n_holes++;
        }

        // Compare recorded and truth traces and count missed surfaces for both
        detray::detail::helix ideal_traj{track, field_view};
        auto [traces_match, n_miss_trace_nav, n_miss_trace_truth, n_error_trace,
              missed_inters] = compare_traces(cfg, truth_trace, recorded_trace,
                                              ideal_traj, i, &(*debug_file));
        // Comparison failed
        if (!traces_match) {
            // Write debug info to file
            *debug_file << "TEST TRACK " << i << ":\n\n"
                        << nav_printer.to_string() << step_printer.to_string();

            // Create SVGs for traces with mismatches
            if (n_miss_trace_truth.n_total() != 0u ||
                n_miss_trace_nav.n_total() != 0u || n_error_trace != 0u) {
                // Only dump SVG if navigator missed a sf. or an error occured
                if (!cfg.display_only_missed() ||
                    (n_miss_trace_nav.n_total() != 0u)) {
                    detray::detector_scanner::display_error(
                        ctx, det, names, cfg.name(), ideal_traj, truth_trace,
                        cfg.svg_style(), i, n_samples, recorded_trace,
                        cfg.verbose());
                }
            }
        }

        // After dummy records insertion, traces should have the same size
        if (truth_trace.size() != recorded_trace.size()) {
            throw std::runtime_error(
                "ERROR: Trace comparison failed: Trace do not have the same "
                "size after intersection matching");
        }

        // Internal error during trace matching
        if (n_error_trace != 0u) {
            throw std::runtime_error(
                "FATAL: Error during track comparison. Please check log "
                "files");
        }

        // Add to global statistics
        n_miss_nav += n_miss_trace_nav;
        n_miss_truth += n_miss_trace_truth;
        n_matching_error += n_error_trace;

        // Count the number of surfaces per type in all traces
        for (std::size_t j = 0; j < truth_trace.size(); ++j) {
            n_surfaces.count(truth_trace[j].intersection.sf_desc,
                             cfg.verbose());
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

        // Maximum number of missed/extra surfaces over all tracks
        trk_stats.n_max_missed_per_trk = math::max(
            trk_stats.n_max_missed_per_trk, n_miss_trace_nav.n_sensitives);
        trk_stats.n_max_extra_per_trk = math::max(trk_stats.n_max_extra_per_trk,
                                                  n_miss_trace_truth.n_total());

        // Count the number of tracks that were tested surccessfully
        trk_stats.n_tracks++;
    }

    // Print results
    const auto n_tracks{static_cast<double>(trk_stats.n_tracks)};
    const auto n_tracks_w_holes{
        static_cast<double>(trk_stats.n_tracks_w_holes)};
    const auto n_tracks_w_extra{
        static_cast<double>(trk_stats.n_tracks_w_extra)};
    const auto n_good_tracks{static_cast<double>(trk_stats.n_good_tracks)};
    const auto n_max_holes_per_trk{
        static_cast<double>(trk_stats.n_max_missed_per_trk)};
    const auto n_max_extra_per_trk{
        static_cast<double>(trk_stats.n_max_extra_per_trk)};

    // Self check
    if (n_holes > n_tracks_w_holes) {
        throw std::runtime_error(
            "Number of tracks with holes is underestimated!");
    }

    // Calculate and display the surface finding efficiency
    if (cfg.verbose()) {
        print_efficiency(trk_stats.n_tracks, n_surfaces, n_miss_nav,
                         n_miss_truth, n_fatal_error, n_matching_error);
    }

    // Column width in output
    constexpr int cw{35};

    std::cout << std::left << std::setw(cw)
              << "No. Tracks with miss. surfaces: " << n_tracks_w_holes << "/"
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
              << "Max no. of miss. surfaces per track: " << n_max_holes_per_trk
              << " (Mean: "
              << ((n_tracks_w_holes == 0)
                      ? "-"
                      : std::to_string(
                            static_cast<double>(n_miss_nav.n_sensitives) /
                            n_tracks_w_holes))
              << ")" << std::endl;
    std::cout << std::left << std::setw(cw + 5)
              << "Max no. of add. surfaces per track: " << n_max_extra_per_trk
              << " (Mean: "
              << ((n_tracks_w_extra == 0)
                      ? "-"
                      : std::to_string(
                            static_cast<double>(n_miss_truth.n_total()) /
                            n_tracks_w_extra))
              << ")" << std::endl;
    std::cout << "\n-----------------------------------" << std::endl;

    return std::tuple{trk_stats,   n_surfaces, n_miss_nav, n_miss_truth,
                      step_traces, mat_traces, mat_records};
}

}  // namespace detray::navigation_validator
