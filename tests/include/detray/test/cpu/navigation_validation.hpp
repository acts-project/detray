/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/test/common/fixture_base.hpp"
#include "detray/test/common/navigation_validation_config.hpp"
#include "detray/test/common/utils/detector_scan_utils.hpp"
#include "detray/test/common/utils/material_validation_utils.hpp"
#include "detray/test/common/utils/navigation_validation_utils.hpp"
#include "detray/tracks/tracks.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray::test {

/// @brief Test class that runs the navigation check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t, template <typename> class scan_type>
class navigation_validation : public test::fixture_base<> {

    using scalar_t = typename detector_t::scalar_type;
    using algebra_t = typename detector_t::algebra_type;
    using free_track_parameters_t = free_track_parameters<algebra_t>;
    using trajectory_type = typename scan_type<algebra_t>::trajectory_type;
    using intersection_trace_t = typename scan_type<
        algebra_t>::template intersection_trace_type<detector_t>;

    /// Switch between rays and helices
    static constexpr auto k_use_rays{
        std::is_same_v<detail::ray<algebra_t>, trajectory_type>};

    public:
    using fixture_type = test::fixture_base<>;
    using config = navigation_validation_config;

    explicit navigation_validation(
        const detector_t &det, const typename detector_t::name_map &names,
        const config &cfg = {},
        const typename detector_t::geometry_context gctx = {})
        : m_cfg{cfg}, m_gctx{gctx}, m_det{det}, m_names{names} {

        if (!m_cfg.whiteboard()) {
            throw std::invalid_argument("No white board was passed to " +
                                        m_cfg.name() + " test");
        }
    }

    /// Run the check
    void TestBody() override {
        using namespace detray;
        using namespace navigation;

        using intersection_t =
            typename intersection_trace_t::value_type::intersection_type;

        // Runge-Kutta stepper
        using hom_bfield_t = bfield::const_field_t;
        using bfield_t =
            std::conditional_t<k_use_rays, navigation_validator::empty_bfield,
                               hom_bfield_t>;
        using rk_stepper_t =
            rk_stepper<typename hom_bfield_t::view_t, algebra_t,
                       unconstrained_step, stepper_rk_policy,
                       stepping::print_inspector>;
        using line_stepper_t =
            line_stepper<algebra_t, unconstrained_step, stepper_default_policy,
                         stepping::print_inspector>;
        using stepper_t =
            std::conditional_t<k_use_rays, line_stepper_t, rk_stepper_t>;

        bfield_t b_field{};
        if constexpr (!k_use_rays) {
            b_field = bfield::create_const_field(m_cfg.B_vector());
        }

        // Use ray or helix
        const std::string det_name{m_det.name(m_names)};
        const std::string truth_data_name{
            k_use_rays ? det_name + "_ray_scan" : det_name + "_helix_scan"};

        /// Collect some statistics
        std::size_t n_tracks{0u}, n_miss{0u}, n_fatal{0u};

        std::cout << "\nINFO: Fetching data from white board..." << std::endl;
        if (!m_cfg.whiteboard()->exists(truth_data_name)) {
            throw std::runtime_error(
                "White board is empty! Please run detector scan first");
        }
        const auto &intersection_traces =
            m_cfg.whiteboard()->template get<std::vector<intersection_trace_t>>(
                truth_data_name);

        std::size_t n_test_tracks{
            std::min(m_cfg.n_tracks(), intersection_traces.size())};
        std::cout << "\nINFO: Running navigation validation on: " << det_name
                  << "...\n"
                  << std::endl;

        std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
        detray::io::file_handle debug_file{"./navigation_validation.txt",
                                           io_mode};

        // Keep a record of track positions and material along the track
        dvector<dvector<navigation::detail::candidate_record<intersection_t>>>
            recorded_traces{};
        dvector<material_validator::material_record<scalar_t>> mat_records{};

        for (const auto &intersection_trace : intersection_traces) {

            if (n_tracks >= m_cfg.n_tracks()) {
                break;
            }

            // Follow the test trajectory with a track and check, if we find
            // the same volumes and distances along the way
            const auto &start = intersection_trace.front();
            const auto &track = start.track_param;
            trajectory_type test_traj = get_parametrized_trajectory(track);

            // Run the propagation
            auto [success, obj_tracer, mat_trace, nav_printer, step_printer] =
                navigation_validator::record_propagation<stepper_t>(
                    m_gctx, m_det, m_cfg.propagation(), track, b_field);

            if (success) {
                // The navigator does not record the initial track position,
                // add it as a dummy record
                obj_tracer.object_trace.insert(
                    obj_tracer.object_trace.begin(),
                    {track.pos(), track.dir(), start.intersection});

                success &= navigation_validator::compare_traces(
                    intersection_trace, obj_tracer.object_trace, test_traj,
                    n_tracks, n_test_tracks);

                if (not success) {
                    // Count mismatches
                    ++n_miss;
                }
            } else {
                // Propagation did not succeed
                ++n_fatal;
            }

            if (not success) {
                // Write debug info to file
                *debug_file << "TEST TRACK " << n_tracks << ":\n\n"
                            << nav_printer.to_string()
                            << step_printer.to_string();

                detector_scanner::display_error(
                    m_gctx, m_det, m_names, m_cfg.name(), test_traj,
                    intersection_trace, m_cfg.svg_style(), n_tracks,
                    n_test_tracks, obj_tracer.object_trace);
            }

            recorded_traces.push_back(std::move(obj_tracer.object_trace));
            mat_records.push_back(mat_trace);

            EXPECT_TRUE(success);

            ++n_tracks;
        }

        // Calculate and display the result
        navigation_validator::print_efficiency(n_tracks, n_miss, n_fatal);

        // Print track positions for plotting
        std::string prefix{k_use_rays ? "ray_" : "helix_"};
        const auto data_path{
            std::filesystem::path{m_cfg.track_param_file()}.parent_path()};
        const auto trk_path{data_path / (prefix + "navigation_track_pos.csv")};
        const auto math_path{data_path / (prefix + "accumulated_material.csv")};

        navigation_validator::write_tracks(trk_path.string(), recorded_traces);
        material_validator::write_material(math_path.string(), mat_records);
    }

    private:
    /// @returns either the helix or ray corresponding to the input track
    /// parameters @param track
    trajectory_type get_parametrized_trajectory(
        const free_track_parameters_t &track) {
        std::unique_ptr<trajectory_type> test_traj{nullptr};
        if constexpr (k_use_rays) {
            test_traj = std::make_unique<trajectory_type>(track);
        } else {
            test_traj =
                std::make_unique<trajectory_type>(track, &(m_cfg.B_vector()));
        }
        return *(test_traj.release());
    }

    /// The configuration of this test
    config m_cfg;
    /// The geometry context to check
    typename detector_t::geometry_context m_gctx{};
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

template <typename detector_t>
using straight_line_navigation =
    navigation_validation<detector_t, detray::ray_scan>;

template <typename detector_t>
using helix_navigation = navigation_validation<detector_t, detray::helix_scan>;

}  // namespace detray::test
