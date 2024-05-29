/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/test/common/fixture_base.hpp"
#include "detray/test/common/material_validation_config.hpp"
#include "detray/test/common/utils/material_validation_utils.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>
#include <stdexcept>
#include <string>

namespace detray::test {

/// Run the material validation on the host
struct run_material_validation {

    inline static const std::string name{"cpu"};

    template <typename detector_t>
    auto operator()(
        vecmem::memory_resource *host_mr, vecmem::memory_resource *,
        const detector_t &det, const propagation::config &cfg,
        const dvector<free_track_parameters<typename detector_t::algebra_type>>
            &tracks) {

        using scalar_t = typename detector_t::scalar_type;

        typename detector_t::geometry_context gctx{};

        dvector<material_validator::material_record<scalar_t>> mat_records{
            host_mr};
        mat_records.reserve(tracks.size());

        for (const auto &[i, track] : detray::views::enumerate(tracks)) {

            auto [success, mat_record] =
                detray::material_validator::record_material(gctx, det, cfg,
                                                            track);
            mat_records.push_back(mat_record);

            if (!success) {
                std::cerr << "ERROR: Propagation failed for track " << i << ": "
                          << "Material record may be incomplete!" << std::endl;
            }
        }

        return mat_records;
    }
};

/// @brief Test class that runs the material validation for a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed outside this class
template <typename detector_t, typename material_validator_t>
class material_validation_impl : public test::fixture_base<> {

    using scalar_t = typename detector_t::scalar_type;
    using algebra_t = typename detector_t::algebra_type;
    using free_track_parameters_t = free_track_parameters<algebra_t>;
    using material_record_t = material_validator::material_record<scalar_t>;

    public:
    using fixture_type = test::fixture_base<>;
    using config = detray::test::material_validation_config;

    explicit material_validation_impl(
        const detector_t &det, const typename detector_t::name_map &names,
        const config &cfg = {},
        const typename detector_t::geometry_context gctx = {})
        : m_cfg{cfg}, m_gctx{gctx}, m_det{det}, m_names{names} {

        if (!m_cfg.whiteboard()) {
            throw std::invalid_argument("No white board was passed to " +
                                        m_cfg.name() + " test");
        }

        // Name of the material scan data collection
        m_scan_data_name = "material_scan_" + m_det.name(m_names);
        m_track_data_name = "material_scan_tracks_" + m_det.name(m_names);

        // Check that data is available in memory
        if (!m_cfg.whiteboard()->exists(m_scan_data_name)) {
            throw std::invalid_argument(
                "Material validation: Could not find scan data on whiteboard."
                "Please run material scan first.");
        }
        if (!m_cfg.whiteboard()->exists(m_track_data_name)) {
            throw std::invalid_argument(
                "Material validation: Could not find track data on whiteboard."
                "Please run material scan first.");
        }
    }

    /// Run the check
    void TestBody() override {
        using namespace detray;

        // Fetch the input data
        const auto &tracks =
            m_cfg.whiteboard()->template get<dvector<free_track_parameters_t>>(
                m_track_data_name);

        const auto &truth_mat_records =
            m_cfg.whiteboard()->template get<dvector<material_record_t>>(
                m_scan_data_name);

        std::cout << "\nINFO: Running material validation on: "
                  << m_det.name(m_names) << "...\n"
                  << std::endl;

        // Run the propagation on device and record the accumulated material
        auto mat_records = material_validator_t{}(
            &m_host_mr, m_cfg.device_mr(), m_det, m_cfg.propagation(), tracks);

        // One material record per track
        ASSERT_EQ(tracks.size(), mat_records.size());

        // Collect some statistics
        std::size_t n_tracks{0u};
        for (std::size_t i = 0u; i < mat_records.size(); ++i) {

            if (n_tracks >= m_cfg.n_tracks()) {
                break;
            }

            const auto &truth_mat = truth_mat_records[i];
            const auto &recorded_mat = mat_records[i];

            EXPECT_NEAR(truth_mat.sX0, recorded_mat.sX0, m_cfg.tol())
                << "Track " << n_tracks << " (X0 / path)";
            EXPECT_NEAR(truth_mat.tX0, recorded_mat.tX0, m_cfg.tol())
                << "Track " << n_tracks << " (X0 / thickness)";
            EXPECT_NEAR(truth_mat.sL0, recorded_mat.sL0, m_cfg.tol())
                << "Track " << n_tracks << " (L0 / path)";
            EXPECT_NEAR(truth_mat.tL0, recorded_mat.tL0, m_cfg.tol())
                << "Track " << n_tracks << " (L0 / thickness)";

            ++n_tracks;
        }

        std::cout << "-----------------------------------\n"
                  << "Tested " << n_tracks << " tracks\n"
                  << "-----------------------------------\n"
                  << std::endl;

        // Print accumulated material per track
        std::filesystem::path mat_path{m_cfg.material_file()};
        auto file_name{material_validator_t::name + "_" +
                       mat_path.stem().string() + "_" + m_det.name(m_names) +
                       mat_path.extension().string()};

        material_validator::write_material(mat_path.replace_filename(file_name),
                                           mat_records);
    }

    private:
    /// The configuration of this test
    config m_cfg;
    /// Vecmem memory resource for the host allocations
    vecmem::host_memory_resource m_host_mr{};
    /// Name of the input data collections
    std::string m_scan_data_name{""};
    std::string m_track_data_name{""};
    /// The geometry context to check
    typename detector_t::geometry_context m_gctx{};
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

template <typename detector_t>
using material_validation =
    material_validation_impl<detector_t, run_material_validation>;

}  // namespace detray::test
