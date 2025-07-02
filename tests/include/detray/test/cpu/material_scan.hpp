/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/tracks/ray.hpp"

// Detray IO include(s)
#include "detray/io/utils/file_handle.hpp"

// Detray test include(s)
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/framework/whiteboard.hpp"
#include "detray/test/validation/detector_scan_utils.hpp"
#include "detray/test/validation/detector_scanner.hpp"
#include "detray/test/validation/material_validation_utils.hpp"

// System include(s)
#include <ios>
#include <iostream>
#include <memory>
#include <string>

namespace detray::test {

/// @brief Test class that runs the material ray scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class material_scan : public test::fixture_base<> {

    using algebra_t = typename detector_t::algebra_type;
    using point2_t = dpoint2D<algebra_t>;
    using scalar_t = dscalar<algebra_t>;
    using ray_t = detail::ray<algebra_t>;
    using track_generator_t = uniform_track_generator<ray_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename track_generator_t::configuration;

        std::string m_name{"material_scan"};
        trk_gen_config_t m_trk_gen_cfg{};
        /// Name of the output file, containing the complete ray material traces
        std::string m_material_file{"material_scan"};
        /// Perform overlaps removal (needed for detector converted from ACTS)
        bool m_overlaps_removal{true};
        /// Tolerance for overlaps
        float m_overlaps_tol{1e-4f * unit<float>::mm};

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
        trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
        const trk_gen_config_t &track_generator() const {
            return m_trk_gen_cfg;
        }
        const std::string &material_file() const { return m_material_file; }
        bool overlaps_removal() const { return m_overlaps_removal; }
        float overlaps_tol() const { return m_overlaps_tol; }
        /// @}

        /// Setters
        /// @{
        config &name(const std::string n) {
            m_name = n;
            return *this;
        }
        config &material_file(const std::string &f) {
            m_material_file = f;
            return *this;
        }
        config &overlaps_removal(const bool o) {
            m_overlaps_removal = o;
            return *this;
        }
        config &overlaps_tol(const scalar_t tol) {
            m_overlaps_tol = static_cast<float>(tol);
            return *this;
        }
        /// @}
    };

    explicit material_scan(
        const detector_t &det, const typename detector_t::name_map &names,
        const config &cfg = {}, std::shared_ptr<test::whiteboard> wb = nullptr,
        const typename detector_t::geometry_context &gctx = {})
        : m_cfg{cfg},
          m_gctx{gctx},
          m_det{det},
          m_names{names},
          m_whiteboard{std::move(wb)} {

        if (!m_whiteboard) {
            throw std::invalid_argument("No white board was passed to " +
                                        m_cfg.name());
        }
    }

    /// Run the ray scan
    void TestBody() override {

        using material_record_t = material_validator::material_record<scalar_t>;

        std::size_t n_tracks{0u};
        auto ray_generator = track_generator_t(m_cfg.track_generator());

        std::cout << "INFO: Running material scan on: " << m_det.name(m_names)
                  << "\n(" << ray_generator.size() << " rays) ...\n"
                  << std::endl;

        // Trace material per ray
        dvector<free_track_parameters<algebra_t>> tracks{};
        tracks.reserve(ray_generator.size());

        dvector<material_record_t> mat_records{};
        mat_records.reserve(ray_generator.size());

        for (const auto ray : ray_generator) {

            // Record all intersections and surfaces along the ray
            auto intersection_record =
                detector_scanner::run<detray::ray_scan>(m_gctx, m_det, ray);

            // Remove certain allowed duplications
            if (m_cfg.overlaps_removal()) {
                detector_scanner::overlaps_removal(intersection_record,
                                                   m_cfg.overlaps_tol());
            }

            if (intersection_record.empty()) {
                std::cout << "ERROR: Intersection trace empty for ray "
                          << n_tracks << "/" << ray_generator.size() << ": "
                          << ray << std::endl;
                break;
            }

            // Record track parameters
            tracks.push_back({ray.pos(), 0.f, ray.dir(), 0.f});

            // New material record
            material_record_t mat_record{};
            mat_record.eta = vector::eta(ray.dir());
            mat_record.phi = vector::phi(ray.dir());

            // Record material for this ray
            for (const auto &[i, record] :
                 detray::views::enumerate(intersection_record)) {

                // Check whether this record has a successor
                if (i < intersection_record.size() - 1) {

                    const auto &current_intr = record.intersection;
                    const auto &next_intr =
                        intersection_record[i + 1].intersection;

                    const bool is_same_intrs{next_intr == current_intr};
                    const bool current_is_pt{current_intr.sf_desc.is_portal()};
                    const bool next_is_pt{next_intr.sf_desc.is_portal()};

                    const bool is_dummy_record{i == 0};

                    // Prevent double counting of material on adjacent portals
                    // (the navigator automatically skips the exit portal)
                    if ((is_same_intrs && next_is_pt && current_is_pt) ||
                        is_dummy_record) {
                        continue;
                    }
                }

                // Don't count the last portal, because navigation terminates
                // before the material is counted
                if (detail::is_invalid_value(record.intersection.volume_link)) {
                    continue;
                }

                const auto sf =
                    geometry::surface{m_det, record.intersection.sf_desc};

                if (!sf.has_material()) {
                    continue;
                }

                const auto &p = record.intersection.local;
                const auto mat_params = sf.template visit_material<
                    material_validator::get_material_params>(
                    point2_t{p[0], p[1]}, cos_angle(m_gctx, sf, ray.dir(), p));

                const scalar_t seg{mat_params.path};
                const scalar_t t{mat_params.thickness};
                const scalar_t mx0{mat_params.mat_X0};
                const scalar_t ml0{mat_params.mat_L0};

                if (mx0 > 0.f) {
                    mat_record.sX0 += seg / mx0;
                    mat_record.tX0 += t / mx0;
                } else {
                    std::cout << "WARNING: Encountered invalid X_0: " << mx0
                              << "\nOn surface: " << sf << std::endl;
                }
                if (ml0 > 0.f) {
                    mat_record.sL0 += seg / ml0;
                    mat_record.tL0 += t / ml0;
                } else {
                    std::cout << "WARNING: Encountered invalid L_0: " << ml0
                              << "\nOn surface: " << sf << std::endl;
                }
            }

            if (mat_record.sX0 == 0.f || mat_record.sL0 == 0.f ||
                mat_record.tX0 == 0.f || mat_record.tL0 == 0.f) {
                std::cout << "WARNING: No material recorded for ray "
                          << n_tracks << "/" << ray_generator.size() << ": "
                          << ray << std::endl;
            }

            mat_records.push_back(mat_record);

            ++n_tracks;
        }

        std::cout << "-----------------------------------\n"
                  << "Tested " << n_tracks << " tracks\n"
                  << "-----------------------------------\n"
                  << std::endl;

        // Write recorded material to csv file
        const std::string det_name{m_det.name(m_names)};
        const auto material_path{std::filesystem::path{m_cfg.material_file()}};
        const auto data_path{material_path.parent_path()};

        std::string file_name{det_name + "_" + material_path.stem().string() +
                              ".csv"};
        material_validator::write_material(data_path / file_name, mat_records);

        // Pin data to whiteboard
        m_whiteboard->add(det_name + "_material_scan", std::move(mat_records));
        m_whiteboard->add(det_name + "_material_scan_tracks",
                          std::move(tracks));
    }

    private:
    /// The configuration of this test
    config m_cfg;
    /// The geometry context to scan
    typename detector_t::geometry_context m_gctx{};
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
    /// Whiteboard to pin data
    std::shared_ptr<test::whiteboard> m_whiteboard{nullptr};
};

}  // namespace detray::test
