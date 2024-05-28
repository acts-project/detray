/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/utils/file_handle.hpp"
#include "detray/materials/detail/material_accessor.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/common/fixture_base.hpp"
#include "detray/test/common/types.hpp"
#include "detray/test/common/utils/detector_scanner.hpp"

// System include(s)
#include <ios>
#include <iostream>
#include <string>

namespace detray {

/// @brief Test class that runs the material ray scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class material_scan : public test::fixture_base<> {

    using algebra_t = typename detector_t::algebra_type;
    using point2_t = typename detector_t::point2_type;
    using scalar_t = typename detector_t::scalar_type;
    using ray_t = detail::ray<algebra_t>;
    using track_generator_t = uniform_track_generator<ray_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename track_generator_t::configuration;

        std::string m_name{"material_scan"};
        trk_gen_config_t m_trk_gen_cfg{};

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
        trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
        const trk_gen_config_t &track_generator() const {
            return m_trk_gen_cfg;
        }
        /// @}

        /// Setters
        /// @{
        config &name(const std::string n) {
            m_name = n;
            return *this;
        }
        /// @}
    };

    explicit material_scan(
        const detector_t &det, const typename detector_t::name_map &names,
        const config &cfg = {},
        const typename detector_t::geometry_context &gctx = {})
        : m_cfg{cfg}, m_gctx{gctx}, m_det{det}, m_names{names} {}

    /// Run the ray scan
    void TestBody() override {

        const typename detector_t::geometry_context gctx{};

        std::size_t n_tracks{0u};
        auto ray_generator = track_generator_t(m_cfg.track_generator());

        // Csv output file
        std::string file_name{m_cfg.name() + "_" + m_det.name(m_names)};
        detray::io::file_handle outfile{
            file_name, ".csv",
            std::ios::out | std::ios::binary | std::ios::trunc};
        *outfile << "eta,phi,mat_sX0,mat_sL0,mat_tX0,mat_tL0" << std::endl;

        std::cout << "INFO: Running material scan on: " << m_det.name(m_names)
                  << "\n(" << ray_generator.size() << " rays) ...\n"
                  << std::endl;

        scalar_t eta{}, phi{}, mat_sX0{}, mat_sL0{}, mat_tX0{}, mat_tL0{};
        for (const auto ray : ray_generator) {

            // Record all intersections and surfaces along the ray
            const auto intersection_record =
                detector_scanner::run<ray_scan>(m_gctx, m_det, ray);

            if (intersection_record.empty()) {
                std::cout << "ERROR: Intersection trace empty for ray "
                          << n_tracks << "/" << ray_generator.size() << ": "
                          << ray << std::endl;
                break;
            }

            eta = getter::eta(ray.dir());
            phi = getter::phi(ray.dir());
            // Total accumulated path length in X0 and L0, respectively
            mat_sX0 = 0.f;
            mat_sL0 = 0.f;
            // Total accumulated material thickness in X0 and L0, respectively
            mat_tX0 = 0.f;
            mat_tL0 = 0.f;

            // Record material for this ray
            for (const auto &record : intersection_record) {

                const auto sf = surface{m_det, record.intersection.sf_desc};

                if (!sf.has_material()) {
                    continue;
                }

                const auto &p = record.intersection.local;
                const auto [seg, t, mx0, ml0] =
                    sf.template visit_material<get_material_params>(
                        point2_t{p[0], p[1]}, sf.cos_angle(gctx, ray.dir(), p));

                if (mx0 > 0.f) {
                    mat_sX0 += seg / mx0;
                    mat_tX0 += t / mx0;
                } else {
                    std::cout << "WARNING: Encountered invalid X_0: " << mx0
                              << "\nOn surface: " << sf << std::endl;
                }
                if (ml0 > 0.f) {
                    mat_sL0 += seg / ml0;
                    mat_tL0 += t / ml0;
                } else {
                    std::cout << "WARNING: Encountered invalid L_0: " << ml0
                              << "\nOn surface: " << sf << std::endl;
                }
            }

            if (mat_sX0 == 0.f or mat_sL0 == 0.f or mat_tX0 == 0.f or
                mat_tL0 == 0.f) {
                std::cout << "WARNING: No material recorded for ray "
                          << n_tracks << "/" << ray_generator.size() << ": "
                          << ray << std::endl;
            }

            *outfile << eta << "," << phi << "," << mat_sX0 << "," << mat_sL0
                     << "," << mat_tX0 << "," << mat_tL0 << std::endl;

            ++n_tracks;
        }
    }

    private:
    /// @brief Functor to retrieve the material parameters of a given
    /// intersection
    struct get_material_params {

        template <typename mat_group_t, typename index_t>
        inline auto operator()(
            [[maybe_unused]] const mat_group_t &mat_group,
            [[maybe_unused]] const index_t &index,
            [[maybe_unused]] const point2_t &loc,
            [[maybe_unused]] const scalar_t cos_inc_angle) const {

            using material_t = typename mat_group_t::value_type;

            constexpr auto inv{detail::invalid_value<scalar_t>()};

            // Access homogeneous surface material or material maps
            if constexpr ((detail::is_hom_material_v<material_t> &&
                           !std::is_same_v<material_t, material<scalar_t>>) ||
                          detail::is_material_map_v<material_t>) {

                // Slab or rod
                const auto mat =
                    detail::material_accessor::get(mat_group, index, loc);

                // Empty material can occur in material maps, skip it
                if (!mat) {
                    // Set the pathlength and thickness to zero so that they
                    // are not counted
                    return std::tuple(scalar_t{0}, scalar_t{0}, inv, inv);
                }

                const scalar_t seg{mat.path_segment(cos_inc_angle, loc[0])};
                const scalar_t t{mat.thickness()};
                const scalar_t mat_X0{mat.get_material().X0()};
                const scalar_t mat_L0{mat.get_material().L0()};

                return std::tuple(seg, t, mat_X0, mat_L0);
            } else {
                return std::tuple(inv, inv, inv, inv);
            }
        }
    };

    /// The configuration of this test
    config m_cfg;
    /// The geometry context to scan
    typename detector_t::geometry_context m_gctx{};
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

}  // namespace detray
