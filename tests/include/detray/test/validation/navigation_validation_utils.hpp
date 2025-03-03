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

// Detray IO include(s)
#include "detray/io/utils/file_handle.hpp"

// Detray test include(s)
#include "detray/test/utils/inspectors.hpp"
#include "detray/test/validation/material_validation_utils.hpp"
#include "detray/test/validation/step_tracer.hpp"

// System include(s)
#include <algorithm>
#include <iomanip>
#include <iterator>
#include <memory>
#include <sstream>

namespace detray::navigation_validator {

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
    template <typename sf_descriptor_t>
    bool count(const sf_descriptor_t &sf_desc) {
        switch (sf_desc.id()) {
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
        };
    }

    /// Count up
    surface_stats &operator+=(const surface_stats &other) {
        n_portals += other.n_portals;
        n_sensitives += other.n_sensitives;
        n_passives += other.n_passives;

        return *this;
    }
};

/// Run the propagation and record test data along the way
template <typename stepper_t, typename detector_t,
          typename bfield_t = empty_bfield>
inline auto record_propagation(
    const typename detector_t::geometry_context ctx,
    vecmem::memory_resource *host_mr, const detector_t &det,
    const propagation::config &cfg,
    const free_track_parameters<typename detector_t::algebra_type> &track,
    const bfield_t &bfield = {}) {

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
    using actor_chain_t =
        actor_chain<pathlimit_aborter_t, step_tracer_t, material_tracer_t>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    // Propagator
    propagator_t prop{cfg};

    // Build actor and propagator states
    typename pathlimit_aborter_t::state pathlimit_aborter_state{
        cfg.stepping.path_limit};
    typename step_tracer_t::state step_tracer_state{*host_mr};
    typename material_tracer_t::state mat_tracer_state{*host_mr};
    auto actor_states = detray::tie(pathlimit_aborter_state, step_tracer_state,
                                    mat_tracer_state);

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

    // Run the propagation
    bool success = prop.propagate(*propagation, actor_states);

    return std::make_tuple(
        success, std::move(obj_tracer),
        std::move(step_tracer_state).release_step_data(),
        std::move(mat_tracer_state).release_material_record(),
        std::move(mat_tracer_state).release_material_steps(),
        std::move(nav_printer), std::move(step_printer));
}

/// Compare the recorded intersection trace to the truth trace
template <typename truth_trace_t, typename recorded_trace_t, typename traj_t>
auto compare_traces(truth_trace_t &truth_trace,
                    recorded_trace_t &recorded_trace, const traj_t &traj,
                    std::size_t trk_no, std::size_t total_n_trks,
                    std::fstream *debug_file = nullptr) {

    using nav_record_t = typename recorded_trace_t::value_type;
    using truth_record_t = typename truth_trace_t::value_type;
    using intersection_t = typename truth_record_t::intersection_type;

    std::stringstream debug_stream;
    std::stringstream matching_stream;
    std::size_t n_inters_nav{recorded_trace.size()};
    std::size_t max_entries{math::max(n_inters_nav, truth_trace.size())};

    // Fill the debug stream with the information from both traces
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

    // Check every single recorded intersection
    surface_stats n_missed_nav{};
    surface_stats n_missed_truth{};
    std::size_t n_errors{0u};
    std::vector<intersection_t> missed_intersections{};
    // Iterate until 'max_entries', because dummy records will be added to the
    // shorter trace
    for (long i = 0; i < static_cast<long>(max_entries); ++i) {

        const auto idx{static_cast<std::size_t>(i)};
        const auto &nav_inters =
            recorded_trace[idx].intersection.sf_desc.barcode();
        const auto &truth_inters =
            truth_trace[idx].intersection.sf_desc.barcode();

        const bool found_same_surfaces{nav_inters == truth_inters};

        if (!found_same_surfaces) {
            // Intersection records at portal boundary might be flipped
            // (the portals overlap completely)
            auto is_swapped_portals = [&recorded_trace,
                                       &truth_trace](const long j) {
                const auto idx_j{static_cast<std::size_t>(j)};

                const auto &current_nav_inters =
                    recorded_trace[idx_j].intersection.sf_desc.barcode();
                const auto &current_truth_inters =
                    truth_trace[idx_j].intersection.sf_desc.barcode();

                const std::size_t next_idx{idx_j + 1u};
                if (next_idx < truth_trace.size() &&
                    next_idx < recorded_trace.size()) {

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

            // Match the barcodes
            auto is_matched_nav = [truth_inters](const nav_record_t &nr) {
                return nr.intersection.sf_desc.barcode() == truth_inters;
            };
            auto is_matched_truth = [nav_inters](const truth_record_t &tr) {
                return tr.intersection.sf_desc.barcode() == nav_inters;
            };

            // Check if the portal order is swapped or the surface appears
            // later in the truth/navigation trace (this means one or
            // multiple surfaces were skipped respectively)
            if (is_swapped_portals(i)) {
                // Have already checked the next record
                ++i;
            } else if (auto last_missed_tr = std::ranges::find_if(
                           std::ranges::begin(truth_trace) + i,
                           std::ranges::end(truth_trace), is_matched_truth);
                       last_missed_tr != std::end(truth_trace)) {

                // The navigator missed a(multiple) surface(s)
                auto first_missed = std::begin(truth_trace) + i;
                const auto n_check{std::distance(first_missed, last_missed_tr)};
                // Number of missed surfaces so far
                long missed_pt{static_cast<long>(n_missed_nav.n_portals)};
                long missed_sn{static_cast<long>(n_missed_nav.n_sensitives)};
                long missed_ps{static_cast<long>(n_missed_nav.n_passives)};

                // Check and record surfaces that were missed: Insert dummy
                // records until traces align again
                for (long j = i; j < i + n_check; ++j) {
                    const auto &truth_sfi =
                        truth_trace[static_cast<std::size_t>(j)].intersection;

                    // Portals may be swapped and wrongfully included in the
                    // range of missed surfaces - skip them
                    if (truth_sfi.sf_desc.is_portal() &&
                        is_swapped_portals(j)) {
                        ++j;
                        continue;
                    }
                    missed_intersections.push_back(truth_sfi);
                    // Insert dummy records to match the truth trace size
                    recorded_trace.insert(recorded_trace.begin() + i,
                                          nav_record_t{});

                    // Count per surface type
                    const bool valid{n_missed_nav.count(truth_sfi.sf_desc)};
                    if (!valid) {
                        matching_stream << "ERROR: Encountered surface of "
                                           "unknown type in intersection "
                                           "trace. Validate the geometry";
                        ++n_errors;
                    }
                }

                auto print_err = [&matching_stream, max_entries, i](
                                     const std::string &sf_type, long n_sf) {
                    matching_stream << "\nERROR: Detray navigator missed "
                                    << n_sf << " " << sf_type << "(s) at: " << i
                                    << "/" << max_entries << ": "
                                    << " (Inserted dummy record(s))";
                };

                // Number of misses on this track
                missed_pt =
                    static_cast<long>(n_missed_nav.n_portals) - missed_pt;
                missed_sn =
                    static_cast<long>(n_missed_nav.n_sensitives) - missed_sn;
                missed_ps =
                    static_cast<long>(n_missed_nav.n_passives) - missed_ps;

                assert(missed_pt >= 0);
                assert(missed_sn >= 0);
                assert(missed_ps >= 0);

                if (missed_pt > 0) {
                    print_err("portal", missed_pt);
                }
                if (missed_sn > 0) {
                    print_err("sensitive", missed_sn);
                }
                if (missed_ps > 0) {
                    print_err("passive", missed_ps);
                }

                // Continue checking where trace might match again
                const long n{missed_pt + missed_sn + missed_ps};
                i += (n - 1);

            } else if (auto last_missed_nav = std::ranges::find_if(
                           std::ranges::begin(recorded_trace) + i,
                           std::ranges::end(recorded_trace), is_matched_nav);
                       last_missed_nav != std::end(recorded_trace)) {
                // The detector scanner missed a(multiple) surface(s)
                auto first_missed = std::begin(recorded_trace) + i;
                const auto n_check{
                    std::distance(first_missed, last_missed_nav)};
                // Number of missed surfaces so far
                long missed_pt{static_cast<long>(n_missed_truth.n_portals)};
                long missed_sn{static_cast<long>(n_missed_truth.n_sensitives)};
                long missed_ps{static_cast<long>(n_missed_truth.n_passives)};

                // Check and record surfaces that were missed
                for (long j = i; j < i + n_check; ++j) {
                    const auto &rec_sfi =
                        recorded_trace[static_cast<std::size_t>(j)]
                            .intersection;

                    // Portals may be swapped and wrongfully included in the
                    // range of missed surfaces - skip them
                    if (rec_sfi.sf_desc.is_portal() && is_swapped_portals(j)) {
                        ++j;
                        continue;
                    }
                    missed_intersections.push_back(rec_sfi);
                    // Insert dummy records to match the truth trace size
                    truth_trace.insert(truth_trace.begin() + i,
                                       truth_record_t{});

                    // Count per surface type
                    const bool valid{n_missed_truth.count(rec_sfi.sf_desc)};
                    if (!valid) {
                        matching_stream << "ERROR: Encountered surface of "
                                           "unknown type in intersection "
                                           "trace. Validate the geometry";
                        ++n_errors;
                    }
                }

                auto print_err = [&matching_stream, max_entries, i](
                                     const std::string &sf_type, long n_sf) {
                    matching_stream << "\nERROR: Detray navigator found "
                                    << n_sf << " additional " << sf_type
                                    << "(s) at: " << i << "/" << max_entries
                                    << " (Inserted dummy record(s))";
                };

                // Number of misses on this track
                missed_pt =
                    static_cast<long>(n_missed_truth.n_portals) - missed_pt;
                missed_sn =
                    static_cast<long>(n_missed_truth.n_sensitives) - missed_sn;
                missed_ps =
                    static_cast<long>(n_missed_truth.n_passives) - missed_ps;

                assert(missed_pt >= 0);
                assert(missed_sn >= 0);
                assert(missed_ps >= 0);

                if (missed_pt > 0) {
                    print_err("portal", missed_pt);
                }
                if (missed_sn > 0) {
                    print_err("sensitive", missed_sn);
                }
                if (missed_ps > 0) {
                    print_err("passive", missed_ps);
                }

                // Continue checking where trace might match again
                const long n{missed_pt + missed_sn + missed_ps};
                i += (n - 1);
            } else {
                // None of the above: Error!
                matching_stream << "\nERROR: More than one consecutive "
                                   "surfaces is unmatched! "
                                << i << "/" << max_entries;

                ++n_errors;
                break;
            }
        }
    }

    const bool any_error{(n_missed_nav.n_total() != 0u) ||
                         (n_missed_truth.n_total() != 0u) || (n_errors != 0u)};

    // Fail the test with some extra information
    EXPECT_TRUE(!any_error)
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
    if (n_missed_nav.n_total() > 1u) {
        std::cout << "WARNING: Detray navigator skipped multiple surfaces: "
                  << n_missed_nav.n_total() << "\n"
                  << std::endl;
    }
    if (n_missed_truth.n_total() > 1u) {
        std::cout << "WARNING: Detray navigator found multiple extra surfaces: "
                  << n_missed_truth.n_total() << "\n"
                  << std::endl;
    }
    // Unknown error occured during matching
    EXPECT_TRUE(n_errors == 0u)
        << "ERROR: Errors during matching: " << n_errors;

    // After inserting the placeholders, do a final check on the trace sizes
    const bool is_size{recorded_trace.size() == truth_trace.size()};
    EXPECT_TRUE(is_size)
        << "ERROR: Intersection traces have different number "
           "of surfaces after matching! Please check unmatched elements\n"
        << debug_stream.str();

    if (!is_size || (n_missed_nav.n_total() != 0u) ||
        (n_missed_truth.n_total() != 0u) || (n_errors != 0u)) {
        return std::make_tuple(false, n_missed_nav, n_missed_truth, n_errors,
                               missed_intersections);
    }

    return std::make_tuple(true, n_missed_nav, n_missed_truth, n_errors,
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
            const auto sf = tracking_surface{det, missed_sfi.sf_desc.barcode()};
            const auto vol = tracking_volume{det, sf.volume()};

            const auto dist = sf.min_dist_to_boundary(missed_sfi.local);
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
/// @NOTE: WIP
inline auto print_efficiency(std::size_t n_tracks,
                             const surface_stats &n_surfaces,
                             const surface_stats &n_miss_nav,
                             const surface_stats &n_miss_truth,
                             std::size_t n_fatal,
                             std::size_t n_matching_error) {
    // Column width in output
    constexpr int cw{20};

    // Print general information
    if (n_miss_nav.n_total() > 0u || n_miss_truth.n_total() > 0u ||
        n_fatal > 0u || n_matching_error > 0u) {

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
            << "\n\nFatal propagation failures:   " << n_fatal
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
    print_eff("Portal sf.", n_surfaces.n_portals, n_miss_nav.n_portals);
    print_eff("Sensitive sf.", n_surfaces.n_sensitives,
              n_miss_nav.n_sensitives);
    print_eff("Passive sf.", n_surfaces.n_passives, n_miss_nav.n_passives);

    std::cout << std::endl;
    print_eff("Surface", n_surfaces.n_total(), n_miss_nav.n_total());

    std::cout << "\n-----------------------------------\n" << std::endl;
}

}  // namespace detray::navigation_validator
