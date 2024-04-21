/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/test/common/fixture_base.hpp"
#include "detray/test/common/navigation_validation_config.hpp"
#include "detray/test/common/utils/detector_scan_utils.hpp"
#include "detray/test/common/utils/detector_scanner.hpp"
#include "detray/test/common/utils/navigation_validation_utils.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/inspectors.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

namespace detray::cuda {

/// Launch the navigation validation kernel
///
/// @param[in] det_view the detector vecmem view
/// @param[in] cfg the propagation configuration
/// @param[in] field_data the magentic field view (maybe an empty field)
/// @param[in] navigation_cache_view the navigation cache vecemem view
/// @param[in] truth_intersection_traces_view vecemem view of the truth data
/// @param[out] recorded_intersections_view vecemem view of the intersections
///                                         recorded by the navigator
template <typename bfield_t, typename detector_t,
          typename intersection_record_t>
void navigation_validation_device(
    typename detector_t::view_type det_view, const propagation::config &cfg,
    bfield_t field_data,
    vecmem::data::jagged_vector_view<
        typename intersection_record_t::intersection_type>
        &navigation_cache_view,
    vecmem::data::jagged_vector_view<const intersection_record_t>
        &truth_intersection_traces_view,
    vecmem::data::jagged_vector_view<navigation::detail::candidate_record<
        typename intersection_record_t::intersection_type>>
        &recorded_intersections_view);

/// Prepare data for device navigation run
template <typename bfield_t, typename detector_t,
          typename intersection_record_t>
inline auto run_navigation_validation(
    vecmem::memory_resource *host_mr, vecmem::memory_resource *dev_mr,
    const detector_t &det, const propagation::config &cfg, bfield_t field_data,
    const std::vector<std::vector<intersection_record_t>>
        &truth_intersection_traces)
    -> vecmem::jagged_vector<navigation::detail::candidate_record<
        typename intersection_record_t::intersection_type>> {

    using intersection_t = typename intersection_record_t::intersection_type;

    // Helper object for performing memory copies (to CUDA devices)
    vecmem::cuda::copy cuda_cpy;

    // Copy the detector to device and get its view
    auto det_buffer = detray::get_buffer(det, *dev_mr, cuda_cpy);
    auto det_view = detray::get_data(det_buffer);

    // Allocate memory for the navigation cache on the device
    const std::size_t n_tracks{truth_intersection_traces.size()};
    auto navigation_cache_buffer =
        detray::create_candidates_buffer(det, n_tracks, *dev_mr, host_mr);
    cuda_cpy.setup(navigation_cache_buffer);

    // Move truth intersection traces data to device
    auto truth_intersection_traces_data =
        vecmem::get_data(truth_intersection_traces, host_mr);
    auto truth_intersection_traces_buffer =
        cuda_cpy.to(truth_intersection_traces_data, *dev_mr, host_mr,
                    vecmem::copy::type::host_to_device);
    vecmem::data::jagged_vector_view<const intersection_record_t>
        truth_intersection_traces_view =
            vecmem::get_data(truth_intersection_traces_buffer);

    // Buffer for the intersections recorded by the navigator
    std::vector<std::size_t> capacities;
    for (const auto &trace : truth_intersection_traces) {
        // Increase the capacity, in case the navigator finds more surfaces
        // than the truth intersections (usually just one)
        capacities.push_back(trace.size() + 10u);
    }

    vecmem::data::jagged_vector_buffer<
        navigation::detail::candidate_record<intersection_t>>
        recorded_intersections_buffer(capacities, *dev_mr, host_mr,
                                      vecmem::data::buffer_type::resizable);
    cuda_cpy.setup(recorded_intersections_buffer);
    auto recorded_intersections_view =
        vecmem::get_data(recorded_intersections_buffer);

    // Run the navigation validation test on device
    navigation_validation_device<bfield_t, detector_t, intersection_record_t>(
        det_view, cfg, field_data, navigation_cache_buffer,
        truth_intersection_traces_view, recorded_intersections_view);

    // Get the result back to the host and pass it on to the checking
    vecmem::jagged_vector<navigation::detail::candidate_record<intersection_t>>
        recorded_intersections(host_mr);
    cuda_cpy(recorded_intersections_buffer, recorded_intersections);

    return recorded_intersections;
}

/// @brief Test class that runs the navigation validation for a given detector
/// on device.
///
/// @note The lifetime of the detector needs to be guaranteed outside this class
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
    using config = detray::test::navigation_validation_config;

    explicit navigation_validation(
        const detector_t &det, const typename detector_t::name_map &names,
        const config &cfg = {},
        const typename detector_t::geometry_context gctx = {})
        : m_cfg{cfg}, m_gctx{gctx}, m_det{det}, m_names{names} {

        if (!m_cfg.whiteboard()) {
            throw std::invalid_argument("No white board was passed to " +
                                        m_cfg.name() + " test");
        }

        // Use ray or helix
        const std::string det_name{m_det.name(m_names)};
        m_truth_data_name = k_use_rays ? det_name + "_ray_scan_for_cuda"
                                       : det_name + "_helix_scan_for_cuda";

        // Pin the data onto the whiteboard
        if (!m_cfg.whiteboard()->exists(m_truth_data_name) &&
            io::file_exists(m_cfg.intersection_file()) &&
            io::file_exists(m_cfg.track_param_file())) {

            // Name clash: Choose alternative name
            if (m_cfg.whiteboard()->exists(m_truth_data_name)) {
                m_truth_data_name = io::alt_file_name(m_truth_data_name);
            }

            std::vector<intersection_trace_t> intersection_traces;

            std::cout << "\nINFO: Reading data from file..." << std::endl;

            // Fill the intersection traces from file
            detray::detector_scanner::read(m_cfg.intersection_file(),
                                           m_cfg.track_param_file(),
                                           intersection_traces);

            m_cfg.whiteboard()->add(m_truth_data_name,
                                    std::move(intersection_traces));
        } else if (m_cfg.whiteboard()->exists(m_truth_data_name)) {
            std::cout << "\nINFO: Fetching data from white board..."
                      << std::endl;
        } else {
            throw std::invalid_argument(
                "Navigation validation: Could not find data files");
        }

        // Check that data is ready
        if (!m_cfg.whiteboard()->exists(m_truth_data_name)) {
            throw std::invalid_argument(
                "Data for navigation check is not on the whiteboard");
        }
    }

    /// Run the check
    void TestBody() override {
        using namespace detray;
        using namespace navigation;

        // Runge-Kutta stepper
        using hom_bfield_t = bfield::const_field_t;
        using bfield_view_t =
            std::conditional_t<k_use_rays, navigation_validator::empty_bfield,
                               hom_bfield_t::view_t>;
        using bfield_t =
            std::conditional_t<k_use_rays, navigation_validator::empty_bfield,
                               hom_bfield_t>;

        bfield_t b_field{};
        if constexpr (!k_use_rays) {
            b_field = bfield::create_const_field(m_cfg.B_vector());
        }

        // Fetch the truth data
        const auto &truth_intersection_traces =
            m_cfg.whiteboard()->template get<std::vector<intersection_trace_t>>(
                m_truth_data_name);

        std::size_t n_test_tracks{
            std::min(m_cfg.n_tracks(), truth_intersection_traces.size())};
        std::cout << "\nINFO: Running device navigation validation on: "
                  << m_det.name(m_names) << "...\n"
                  << std::endl;

        // Run the propagation on device and record the navigation data
        auto recorded_intersections = run_navigation_validation<bfield_view_t>(
            &m_host_mr, &m_dev_mr, m_det, m_cfg.propagation(), b_field,
            truth_intersection_traces);

        // Collect some statistics
        std::size_t n_tracks{0u}, n_miss{0u}, n_fatal{0u};

        EXPECT_EQ(recorded_intersections.size(),
                  truth_intersection_traces.size());

        for (std::size_t i = 0u; i < truth_intersection_traces.size(); ++i) {
            const auto &truth_trace = truth_intersection_traces[i];
            const auto &recorded_trace = recorded_intersections[i];

            if (n_tracks >= m_cfg.n_tracks()) {
                break;
            }

            // Get the original test trajectory (ray or helix)
            const auto &trck_param = truth_trace.front().track_param;
            trajectory_type test_traj = get_parametrized_trajectory(trck_param);

            // Recorded only the start position, which added by default
            bool success{true};
            if (truth_trace.size() == 1) {
                // Propagation did not succeed
                success = false;
                ++n_fatal;
            } else {
                // Compare truth and recorded data elementwise
                success &= navigation_validator::compare_traces(
                    truth_trace, recorded_trace, test_traj, n_tracks,
                    n_test_tracks);

                if (not success) {
                    // Count mismatches
                    ++n_miss;
                }
            }

            if (not success) {
                detector_scanner::display_error(
                    m_gctx, m_det, m_names, m_cfg.name(), test_traj,
                    truth_trace, m_cfg.svg_style(), n_tracks, n_test_tracks,
                    recorded_trace);
            }

            EXPECT_TRUE(success);

            ++n_tracks;
        }

        // Calculate and display the result
        navigation_validator::print_efficiency(n_tracks, n_miss, n_fatal);

        // Print track positions for plotting
        std::string prefix{k_use_rays ? "ray_" : "helix_"};
        const auto data_path{
            std::filesystem::path{m_cfg.track_param_file()}.parent_path() /
            (prefix + "navigation_track_pos_cuda.csv")};

        navigation_validator::write_tracks(data_path.string(),
                                           recorded_intersections);
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

    /// Vecmem memory resource for the host allocations
    vecmem::host_memory_resource m_host_mr{};
    /// Vecmem memory resource for the device allocations
    vecmem::cuda::device_memory_resource m_dev_mr{};
    /// The configuration of this test
    config m_cfg;
    /// Name of the truth data collection
    std::string m_truth_data_name{""};
    /// The geometry context to check
    typename detector_t::geometry_context m_gctx{};
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

template <typename detector_t>
using straight_line_navigation =
    detray::cuda::navigation_validation<detector_t, detray::ray_scan>;

template <typename detector_t>
using helix_navigation =
    detray::cuda::navigation_validation<detector_t, detray::helix_scan>;

}  // namespace detray::cuda
